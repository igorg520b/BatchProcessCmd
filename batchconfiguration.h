#ifndef BATCHCONFIGURATION_H
#define BATCHCONFIGURATION_H

#include <QString>
#include <QSettings>
#include <QDebug>
#include <QList>
#include <QStringList>
#include <QVariant>
#include <QDir>

#include <string>

#include "mesh.h"


class BatchConfiguration
{
public:
    BatchConfiguration();


    QList<double> czsStrengths, rotationAngles;
    QStringList mshFileNames;
    QString pathToAbaqus;
    int numberOfCores; // for job

    QString batchFileName;          // where this configuration is saved
    QString BatchName() const;      // same as batchFileName, without extension and path

    void Load(QString fileName);

    static void ParseValueList(QString CSV, QList<double> &outVec);
    static QString ListToString(QList<double> &vec);
    static QString ToShortName(QString filePath);

    // GENERATION OF PY FILES
    void PrepareTable();        // list of .PY files to be generated
    void ProducePYFiles();      // save .PY files into a fodler
    void Convert_PY_to_INP();   // invoke Abaqus

    struct TableEntry
    {
        int id;
        QString mshFileName;
        double czStrength;
        double rotationAngle;
        QString pyFileName;
    };

    QList<TableEntry> tableEntries;

};

#endif // BATCHCONFIGURATION_H
