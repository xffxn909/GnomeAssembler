#include "gnome_assembler/main_window.h"

#include <QApplication>
#include <QPushButton>

#include <iostream>

int main(int argc, char* argv[]) {
    QApplication application(argc, argv);
    gnome::MainWindow window;

    if (window.windowTitle() != "GnomeAssembler - FASTQ Cleaning") {
        std::cerr << "The GUI window title is incorrect.\n";
        return 1;
    }

    bool has_run_button = false;
    for (const QPushButton* button : window.findChildren<QPushButton*>()) {
        if (button->text() == "Clean FASTQ") {
            has_run_button = true;
        }
    }
    if (!has_run_button) {
        std::cerr << "The GUI does not contain its cleaning action button.\n";
        return 1;
    }

    std::cout << "Qt cleaning window constructed successfully.\n";
    return 0;
}
