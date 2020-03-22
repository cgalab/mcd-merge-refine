#include <iostream>
#include <chrono>
#include <boost/program_options.hpp>
#include "decomp.h"
#include "view.h"

namespace po = boost::program_options;

int main(int argc, char *argv[]) {
    po::options_description desc("Allowed options");
    desc.add_options()
            ("help", "produce help message")
            ("input", po::value<std::string>(), "input file")
            ("output", po::value<std::string>(), "output file")
            ("status-fd", po::value<int>()->default_value(-1), "status-fd")
            ("view", "open OpenGL interface")
            ;

    po::positional_options_description p;
    p.add("input", 1);
    p.add("output", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << desc << "\n";
        return 1;
    }

    if (vm.count("input")) {
        const auto inFile = vm["input"].as<std::string>();
        
        std::string outFilePath = "";
        if (vm.count("output")) {
            outFilePath = vm["output"].as<std::string>();
        }
        
        const auto cd = ConvexDecomp(inFile, outFilePath, vm["status-fd"].as<int>());

#ifdef ENABLE_VIEW
        bool bEnableView = vm.count("view");
        
        if (bEnableView) {
            QApplication app(argc, argv);
            QCoreApplication::setApplicationName("Convex");
            QCoreApplication::setApplicationVersion(QT_VERSION_STR);
            MainWindow mainWindow;
            mainWindow.addItem(new ConvexDecompGraphicsItem(cd));
            mainWindow.show();

            return app.exec();
        }
#endif
    }
    
    return 0;
}
