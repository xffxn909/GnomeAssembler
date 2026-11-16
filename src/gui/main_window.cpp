#include "gnome_assembler/main_window.h"

#include "gnome_assembler/cleaning_options.h"
#include "gnome_assembler/output_writer.h"

#include <QDir>
#include <QFileDialog>
#include <QFileInfo>
#include <QFormLayout>
#include <QHBoxLayout>
#include <QLabel>
#include <QLineEdit>
#include <QMessageBox>
#include <QPlainTextEdit>
#include <QPushButton>
#include <QSpinBox>
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

MainWindow::MainWindow(QWidget* parent) : QMainWindow(parent) {
    setWindowTitle("GnomeAssembler - FASTQ Cleaning");
    resize(720, 570);

    auto* content = new QWidget(this);
    auto* layout = new QVBoxLayout(content);
    auto* introduction = new QLabel(
        "Milestone 1: clean Phred+33 FASTQ reads before genome assembly. "
        "Choose a file, set the filters, then run cleaning.", content);
    introduction->setWordWrap(true);
    layout->addWidget(introduction);

    auto* form = new QFormLayout;
    QWidget* reads_row = path_row(reads_path_, "Browse...", content);
    QWidget* output_row = path_row(output_prefix_, "Folder...", content);
    form->addRow("FASTQ input:", reads_row);
    form->addRow("Output prefix:", output_row);
    output_prefix_->setPlaceholderText("Choose a folder or enter a file prefix");
    output_prefix_->setToolTip(
        "An existing folder is accepted. Output files will be named cleaned_reads.cleaned.* inside it.");

    minimum_quality_ = new QSpinBox(content);
    minimum_quality_->setRange(0, 93);
    minimum_quality_->setValue(20);
    minimum_quality_->setToolTip("Trim low-quality bases from the 3' end.");
    form->addRow("Minimum quality:", minimum_quality_);

    minimum_length_ = new QSpinBox(content);
    minimum_length_->setRange(0, 1000000);
    minimum_length_->setValue(50);
    form->addRow("Minimum read length:", minimum_length_);

    adapter_prefix_ = new QLineEdit(content);
    adapter_prefix_->setPlaceholderText("Optional DNA sequence at the 5' end");
    form->addRow("5' adapter:", adapter_prefix_);
    adapter_suffix_ = new QLineEdit(content);
    adapter_suffix_->setPlaceholderText("Optional DNA sequence at the 3' end");
    form->addRow("3' adapter:", adapter_suffix_);
    layout->addLayout(form);

    auto* run_button = new QPushButton("Clean FASTQ", content);
    run_button->setDefault(true);
    layout->addWidget(run_button);
    layout->addWidget(new QLabel("Run log", content));
    log_ = new QPlainTextEdit(content);
    log_->setReadOnly(true);
    layout->addWidget(log_, 1);
    setCentralWidget(content);

    connect(reads_row->findChild<QPushButton*>(), &QPushButton::clicked,
            this, [this] { choose_reads_file(); });
    connect(output_row->findChild<QPushButton*>(), &QPushButton::clicked,
            this, [this] { choose_output_folder(); });
    connect(run_button, &QPushButton::clicked, this, [this] { run_cleaning(); });
}

void MainWindow::choose_reads_file() {
    const QString path = QFileDialog::getOpenFileName(this, "Select FASTQ file", {},
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

void MainWindow::run_cleaning() {
    CleaningOptions options;
    options.reads_path = reads_path_->text().toStdString();
    options.output_prefix = output_prefix_->text().toStdString();
    options.minimum_quality = minimum_quality_->value();
    options.minimum_read_length = static_cast<std::size_t>(minimum_length_->value());
    options.adapter_prefix = adapter_prefix_->text().toStdString();
    options.adapter_suffix = adapter_suffix_->text().toStdString();

    const PipelineResult result = run_cleaning_to_files(options);
    log_->setPlainText(QString::fromStdString(result.summary));
    if (result.error.empty()) {
        log_->appendPlainText("\nFinished successfully.");
    } else {
        log_->appendPlainText("\nError: " + QString::fromStdString(result.error));
        QMessageBox::warning(this, "Cleaning could not finish", QString::fromStdString(result.error));
    }
}

} // namespace gnome
