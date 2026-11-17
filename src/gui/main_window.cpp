#include "gnome_assembler/main_window.h"

#include "gnome_assembler/assembly_pipeline.h"
#include "gnome_assembler/cleaning_options.h"
#include "gnome_assembler/graph_pipeline.h"
#include "gnome_assembler/kmer_options.h"
#include "gnome_assembler/output_writer.h"

#include <QDir>
#include <QFileDialog>
#include <QFileInfo>
#include <QFormLayout>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QLabel>
#include <QLineEdit>
#include <QMessageBox>
#include <QMetaObject>
#include <QPlainTextEdit>
#include <QProgressBar>
#include <QPushButton>
#include <QSpinBox>
#include <QThread>
#include <QVBoxLayout>
#include <QWidget>

namespace gnome {

namespace {

QWidget* path_row(QLineEdit*& field, const QString& button_text, QWidget* parent) {
    auto* row = new QWidget(parent);
    auto* layout = new QHBoxLayout(row);
    layout->setContentsMargins(0, 0, 0, 0);
    field = new QLineEdit(row);
    auto* button = new QPushButton(button_text, row);
    layout->addWidget(field, 1);
    layout->addWidget(button);
    return row;
}

} // namespace

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------

MainWindow::MainWindow(QWidget* parent) : QMainWindow(parent) {
    setWindowTitle("GnomeAssembler - FASTQ Cleaning, Graph & Assembly");
    resize(760, 700);

    auto* content = new QWidget(this);
    auto* layout = new QVBoxLayout(content);

    auto* introduction = new QLabel(
        "GnomeAssembler: clean Phred+33 FASTQ reads and build a de Bruijn graph. "
        "Choose a file, adjust the filters, then clean reads, build a graph, or assemble contigs.", content);
    introduction->setWordWrap(true);
    layout->addWidget(introduction);

    // ---- Input / Output ----
    auto* io_box = new QGroupBox("Input / Output", content);
    auto* io_form = new QFormLayout(io_box);
    QWidget* reads_row = path_row(reads_path_, "Browse...", io_box);
    QWidget* output_row = path_row(output_prefix_, "Folder...", io_box);
    io_form->addRow("FASTQ input:", reads_row);
    io_form->addRow("Output prefix:", output_row);
    output_prefix_->setPlaceholderText("Choose a folder or enter a file prefix");
    output_prefix_->setToolTip(
        "An existing folder is accepted. Output files will be named inside it.");
    layout->addWidget(io_box);

    // ---- Cleaning options ----
    auto* clean_box = new QGroupBox("Cleaning Options", content);
    auto* clean_form = new QFormLayout(clean_box);

    minimum_quality_ = new QSpinBox(clean_box);
    minimum_quality_->setRange(0, 93);
    minimum_quality_->setValue(20);
    minimum_quality_->setToolTip("Trim low-quality bases from the 3' end.");
    clean_form->addRow("Minimum quality:", minimum_quality_);

    minimum_length_ = new QSpinBox(clean_box);
    minimum_length_->setRange(0, 1000000);
    minimum_length_->setValue(50);
    clean_form->addRow("Minimum read length:", minimum_length_);

    adapter_prefix_ = new QLineEdit(clean_box);
    adapter_prefix_->setPlaceholderText("Optional DNA sequence at the 5' end");
    clean_form->addRow("5' adapter:", adapter_prefix_);

    adapter_suffix_ = new QLineEdit(clean_box);
    adapter_suffix_->setPlaceholderText("Optional DNA sequence at the 3' end");
    clean_form->addRow("3' adapter:", adapter_suffix_);

    layout->addWidget(clean_box);

    // ---- K-mer / Graph options ----
    auto* kmer_box = new QGroupBox("K-mer / Graph Options (Milestone 2)", content);
    auto* kmer_form = new QFormLayout(kmer_box);

    kmer_length_ = new QSpinBox(kmer_box);
    kmer_length_->setRange(2, 999);
    kmer_length_->setValue(21);
    kmer_length_->setToolTip("Length k of each k-mer. Use an odd number (e.g. 21, 31).");
    kmer_form->addRow("k-mer length:", kmer_length_);

    min_kmer_abundance_ = new QSpinBox(kmer_box);
    min_kmer_abundance_->setRange(1, 1000000);
    min_kmer_abundance_->setValue(2);
    min_kmer_abundance_->setToolTip(
        "Discard k-mers seen fewer than this many times (removes sequencing errors).");
    kmer_form->addRow("Min k-mer abundance:", min_kmer_abundance_);

    layout->addWidget(kmer_box);

    // ---- Buttons ----
    auto* button_row = new QWidget(content);
    auto* button_layout = new QHBoxLayout(button_row);
    button_layout->setContentsMargins(0, 0, 0, 0);

    clean_button_ = new QPushButton("Clean FASTQ", button_row);
    clean_button_->setDefault(true);
    button_layout->addWidget(clean_button_);

    graph_button_ = new QPushButton("Build Graph", button_row);
    graph_button_->setToolTip("Run cleaning then k-mer counting and graph construction.");
    button_layout->addWidget(graph_button_);

    assemble_button_ = new QPushButton("Assemble Contigs", button_row);
    assemble_button_->setToolTip("Run the shared pipeline and traverse a valid Eulerian graph into contigs.");
    button_layout->addWidget(assemble_button_);

    layout->addWidget(button_row);

    // ---- Progress bar ----
    progress_bar_ = new QProgressBar(content);
    progress_bar_->setRange(0, 0); // indeterminate (marquee)
    progress_bar_->setVisible(false);
    layout->addWidget(progress_bar_);

    // ---- Log ----
    layout->addWidget(new QLabel("Run log", content));
    log_ = new QPlainTextEdit(content);
    log_->setReadOnly(true);
    layout->addWidget(log_, 1);

    setCentralWidget(content);

    // ---- Connections ----
    connect(reads_row->findChild<QPushButton*>(), &QPushButton::clicked,
            this, [this] { choose_reads_file(); });
    connect(output_row->findChild<QPushButton*>(), &QPushButton::clicked,
            this, [this] { choose_output_folder(); });
    connect(clean_button_, &QPushButton::clicked,
            this, [this] { run_cleaning(); });
    connect(graph_button_, &QPushButton::clicked,
            this, [this] { run_graph_pipeline(); });
    connect(assemble_button_, &QPushButton::clicked,
            this, [this] { run_assembly_pipeline(); });
}

MainWindow::~MainWindow() {
    if (worker_thread_ && worker_thread_->isRunning()) {
        worker_thread_->quit();
        worker_thread_->wait();
    }
}

// ---------------------------------------------------------------------------
// File picker helpers
// ---------------------------------------------------------------------------

void MainWindow::choose_reads_file() {
    const QString path = QFileDialog::getOpenFileName(
        this, "Select FASTQ file", {},
        "FASTQ files (*.fastq *.fq *.txt);;All files (*)");
    if (!path.isEmpty()) {
        reads_path_->setText(path);
        if (output_prefix_->text().isEmpty()) {
            const QFileInfo input_file(path);
            output_prefix_->setText(
                QDir(input_file.absolutePath()).filePath("GnomeAssemblerResults/cleaned_reads"));
        }
    }
}

void MainWindow::choose_output_folder() {
    const QString folder = QFileDialog::getExistingDirectory(this, "Select output folder");
    if (!folder.isEmpty()) output_prefix_->setText(folder);
}

// ---------------------------------------------------------------------------
// Worker thread management
// ---------------------------------------------------------------------------

void MainWindow::set_running(bool running) {
    clean_button_->setEnabled(!running);
    graph_button_->setEnabled(!running);
    assemble_button_->setEnabled(!running);
    progress_bar_->setVisible(running);
}

void MainWindow::pipeline_finished(const QString& summary, const QString& error) {
    set_running(false);
    log_->setPlainText(summary);
    if (error.isEmpty()) {
        log_->appendPlainText("\nFinished successfully.");
    } else {
        log_->appendPlainText("\nError: " + error);
        QMessageBox::warning(this, "Pipeline could not finish", error);
    }

    if (worker_thread_) {
        worker_thread_->quit();
        worker_thread_->wait();
        worker_thread_->deleteLater();
        worker_thread_ = nullptr;
    }
}

// ---------------------------------------------------------------------------
// Run cleaning (off-thread)
// ---------------------------------------------------------------------------

void MainWindow::run_cleaning() {
    CleaningOptions options;
    options.reads_path           = reads_path_->text().toStdString();
    options.output_prefix        = output_prefix_->text().toStdString();
    options.minimum_quality      = minimum_quality_->value();
    options.minimum_read_length  = static_cast<std::size_t>(minimum_length_->value());
    options.adapter_prefix       = adapter_prefix_->text().toStdString();
    options.adapter_suffix       = adapter_suffix_->text().toStdString();

    if (options.reads_path.empty() || options.output_prefix.empty()) {
        QMessageBox::warning(this, "Missing input",
                             "Please select a FASTQ file and an output prefix.");
        return;
    }

    log_->clear();
    log_->appendPlainText("Cleaning FASTQ reads\u2026");
    set_running(true);

    worker_thread_ = new QThread(this);
    // Capture options by value; post result back to the GUI thread via queued connection.
    QThread* thread = worker_thread_;
    MainWindow* self = this;
    connect(thread, &QThread::started, [options, self, thread]() {
        const PipelineResult result = run_cleaning_to_files(options);
        const QString summary = QString::fromStdString(result.summary);
        const QString error   = QString::fromStdString(result.error);
        QMetaObject::invokeMethod(self, [self, summary, error]() {
            self->pipeline_finished(summary, error);
        }, Qt::QueuedConnection);
        thread->quit();
    });
    worker_thread_->start();
}

// ---------------------------------------------------------------------------
// Run graph pipeline (off-thread)
// ---------------------------------------------------------------------------

void MainWindow::run_graph_pipeline() {
    GraphPipelineOptions options;
    options.cleaning.reads_path          = reads_path_->text().toStdString();
    options.cleaning.output_prefix       = output_prefix_->text().toStdString();
    options.cleaning.minimum_quality     = minimum_quality_->value();
    options.cleaning.minimum_read_length = static_cast<std::size_t>(minimum_length_->value());
    options.cleaning.adapter_prefix      = adapter_prefix_->text().toStdString();
    options.cleaning.adapter_suffix      = adapter_suffix_->text().toStdString();
    options.kmer.kmer_length             = kmer_length_->value();
    options.kmer.min_kmer_abundance      = static_cast<uint32_t>(min_kmer_abundance_->value());
    options.output_prefix                = options.cleaning.output_prefix;

    if (options.cleaning.reads_path.empty()) {
        QMessageBox::warning(this, "Missing input", "Please select a FASTQ file.");
        return;
    }

    log_->clear();
    log_->appendPlainText("Building de Bruijn graph\u2026");
    set_running(true);

    worker_thread_ = new QThread(this);
    QThread* thread = worker_thread_;
    MainWindow* self = this;
    connect(thread, &QThread::started, [options, self, thread]() {
        const GraphPipelineResult result = gnome::run_graph_pipeline(options);
        const QString summary = QString::fromStdString(result.summary);
        const QString error   = QString::fromStdString(result.error);
        QMetaObject::invokeMethod(self, [self, summary, error]() {
            self->pipeline_finished(summary, error);
        }, Qt::QueuedConnection);
        thread->quit();
    });
    worker_thread_->start();
}

// ---------------------------------------------------------------------------
// Run assembly pipeline (off-thread)
// ---------------------------------------------------------------------------

void MainWindow::run_assembly_pipeline() {
    GraphPipelineOptions options;
    options.cleaning.reads_path          = reads_path_->text().toStdString();
    options.cleaning.output_prefix       = output_prefix_->text().toStdString();
    options.cleaning.minimum_quality     = minimum_quality_->value();
    options.cleaning.minimum_read_length = static_cast<std::size_t>(minimum_length_->value());
    options.cleaning.adapter_prefix      = adapter_prefix_->text().toStdString();
    options.cleaning.adapter_suffix      = adapter_suffix_->text().toStdString();
    options.kmer.kmer_length             = kmer_length_->value();
    options.kmer.min_kmer_abundance      = static_cast<uint32_t>(min_kmer_abundance_->value());
    options.output_prefix                = options.cleaning.output_prefix;

    if (options.cleaning.reads_path.empty()) {
        QMessageBox::warning(this, "Missing input", "Please select a FASTQ file.");
        return;
    }

    log_->clear();
    log_->appendPlainText("Assembling contigs with Eulerian traversal…");
    set_running(true);

    worker_thread_ = new QThread(this);
    QThread* thread = worker_thread_;
    MainWindow* self = this;
    connect(thread, &QThread::started, [options, self, thread]() {
        const AssemblyPipelineResult result = gnome::run_assembly_pipeline(options);
        const QString summary = QString::fromStdString(result.summary);
        const QString error   = QString::fromStdString(result.error);
        QMetaObject::invokeMethod(self, [self, summary, error]() {
            self->pipeline_finished(summary, error);
        }, Qt::QueuedConnection);
        thread->quit();
    });
    worker_thread_->start();
}

} // namespace gnome
