#pragma once

#include <QMainWindow>

class QLineEdit;
class QPlainTextEdit;
class QSpinBox;

namespace gnome {

class MainWindow : public QMainWindow {
public:
    explicit MainWindow(QWidget* parent = nullptr);

private:
    void choose_reads_file();
    void choose_output_folder();
    void run_cleaning();

    QLineEdit* reads_path_ = nullptr;
    QLineEdit* output_prefix_ = nullptr;
    QLineEdit* adapter_prefix_ = nullptr;
    QLineEdit* adapter_suffix_ = nullptr;
    QSpinBox* minimum_quality_ = nullptr;
    QSpinBox* minimum_length_ = nullptr;
    QPlainTextEdit* log_ = nullptr;
};

} // namespace gnome
