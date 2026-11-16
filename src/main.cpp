#include "gnome_assembler/command_line.h"
#include "gnome_assembler/main_window.h"
#include "gnome_assembler/output_writer.h"

#include <QApplication>
#include <iostream>

int main(int argc, char* argv[]) {
    const auto command = gnome::parse_command_line(argc, argv);

    if (command.show_help) {
        std::cout << gnome::command_help();
        return 0;
    }
    if (!command.error.empty()) {
        std::cerr << "Error: " << command.error << "\n\n" << gnome::command_help();
        return 2;
    }
    if (command.has_cleaning_request) {
        const auto result = gnome::run_cleaning_to_files(command.options);
        std::cout << result.summary;
        if (!result.error.empty()) {
            std::cerr << "Error: " << result.error << '\n';
            return 1;
        }
        return 0;
    }

    QApplication application(argc, argv);
    gnome::MainWindow window;
    window.show();
    return application.exec();
}
