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

    QCommandLineOption doNotCreateIndenterOption(QStringList() << "n" << "noindenter",
                QCoreApplication::translate("main", "Do not create the indenter"));
    parser.addOption(doNotCreateIndenterOption);

    QCommandLineOption rhitaOption(QStringList() << "r" << "rhita",
                QCoreApplication::translate("main", "Set up simulation for RHITA"));
    parser.addOption(rhitaOption);


    // Process the actual command line arguments given by the user
    parser.process(a);

    const QStringList args = parser.positionalArguments();
    QString configFileName = args.at(0);
    qDebug() << configFileName;

    bool doNotCreateIndenter = parser.isSet(doNotCreateIndenterOption);
    bool rhitaSetup = parser.isSet(rhitaOption);

    qDebug() << "doNotCreateIndenter" << doNotCreateIndenter;
    qDebug() << "rhitaSetup" << rhitaSetup;

    QFileInfo fi(configFileName);

    QString baseFileName = fi.baseName();
    qDebug() << baseFileName;

    BatchConfiguration bc;
    bc.Load(configFileName);
    bc.PrepareTable();

    gmsh::initialize();
    bc.ProducePYFiles(doNotCreateIndenter, rhitaSetup);
    bc.Convert_PY_to_INP();

//    return a.exec();
}
