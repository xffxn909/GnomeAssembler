#include "gnome_assembler/main_window.h"

#include <QApplication>
#include <QPushButton>
#include <QSpinBox>

#include <iostream>

int main(int argc, char* argv[]) {
    QApplication application(argc, argv);
    gnome::MainWindow window;

    if (window.windowTitle() != "GnomeAssembler - FASTQ Cleaning, Graph & Assembly") {
        std::cerr << "The GUI window title is incorrect: "
                  << window.windowTitle().toStdString() << "\n";
        return 1;
    }

    bool has_clean_button = false;
    bool has_graph_button = false;
    bool has_assemble_button = false;
    bool has_open_output_button = false;
    for (const QPushButton* button : window.findChildren<QPushButton*>()) {
        if (button->text() == "Clean FASTQ") has_clean_button = true;
        if (button->text() == "Build Graph") has_graph_button = true;
        if (button->text() == "Assemble Contigs") has_assemble_button = true;
        if (button->text() == "Open Output Folder") has_open_output_button = true;
    }
    if (!has_clean_button) {
        std::cerr << "The GUI does not contain the Clean FASTQ button.\n";
        return 1;
    }
    if (!has_graph_button) {
        std::cerr << "The GUI does not contain the Build Graph button.\n";
        return 1;
    }
    if (!has_assemble_button) {
        std::cerr << "The GUI does not contain the Assemble Contigs button.\n";
        return 1;
    }
    if (!has_open_output_button) {
        std::cerr << "The GUI does not contain the Open Output Folder button.\n";
        return 1;
    }
    if (window.findChildren<QSpinBox*>().size() < 6) {
        std::cerr << "The GUI does not expose the expected BLAST threshold controls.\n";
        return 1;
    }

    std::cout << "Qt cleaning/graph window constructed successfully.\n";
    return 0;
}
