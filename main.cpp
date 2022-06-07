#include <QCoreApplication>
#include <QCommandLineParser>
#include <QDebug>
#include <QFileInfo>

#include "batchconfiguration.h"
#include "gmsh.h"


int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    QCoreApplication::setApplicationName("BatchProcessorCmd");
    QCoreApplication::setApplicationVersion("1.0");

    QCommandLineParser parser;
    parser.setApplicationDescription("Batch generator of Abaqus .inp files");
    parser.addHelpOption();
    parser.addVersionOption();
    parser.addPositionalArgument("source", QCoreApplication::translate("main", "Configuration file."));

    QCommandLineOption noinpOption(QStringList() << "n" << "noinp",
                QCoreApplication::translate("main", "Don't create the input file"));
    parser.addOption(noinpOption);

    // Process the actual command line arguments given by the user
    parser.process(a);

    bool doNotCreateInputFile = parser.isSet(noinpOption);

    const QStringList args = parser.positionalArguments();
    QString configFileName = args.at(0);
    qDebug() << configFileName;

    QFileInfo fi(configFileName);

    QString baseFileName = fi.baseName();
    qDebug() << baseFileName;

    BatchConfiguration bc;
    bc.Load(configFileName);
    bc.PrepareTable();

    gmsh::initialize();
    bc.ProducePYFiles();

    if(!doNotCreateInputFile) bc.Convert_PY_to_INP();

//    return a.exec();
}
