#pragma once

#include <QMainWindow>

class QLineEdit;
class QPlainTextEdit;
class QPushButton;
class QProgressBar;
class QSpinBox;
class QThread;

namespace gnome {

class MainWindow : public QMainWindow {
    Q_OBJECT
public:
    explicit MainWindow(QWidget* parent = nullptr);
    ~MainWindow() override;

private:
    void choose_reads_file();
    void choose_output_folder();
    void open_output_folder();
    void run_cleaning();
    void run_graph_pipeline();
    void run_assembly_pipeline();
    void set_running(bool running);
    void pipeline_finished(const QString& summary, const QString& error);

    // --- Shared inputs ---
    QLineEdit* reads_path_ = nullptr;
    QLineEdit* output_prefix_ = nullptr;
    QLineEdit* adapter_prefix_ = nullptr;
    QLineEdit* adapter_suffix_ = nullptr;
    QSpinBox* minimum_quality_ = nullptr;
    QSpinBox* minimum_length_ = nullptr;

    // --- Milestone 2 k-mer inputs ---
    QSpinBox* kmer_length_ = nullptr;
    QSpinBox* min_kmer_abundance_ = nullptr;

    // --- Action controls ---
    QPushButton* clean_button_ = nullptr;
    QPushButton* graph_button_ = nullptr;
    QPushButton* assemble_button_ = nullptr;
    QPushButton* open_output_button_ = nullptr;
    QProgressBar* progress_bar_ = nullptr;
    QPlainTextEdit* log_ = nullptr;

    // --- Background thread (owned, nullable) ---
    QThread* worker_thread_ = nullptr;
};

} // namespace gnome
