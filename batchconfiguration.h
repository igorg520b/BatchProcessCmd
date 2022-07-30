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
    BatchConfiguration() = default;


    QString mshFileName;
    int dim = 3;


    double czsStrength = 4e5;
    double YoungsModulus = 9e9;
    double czElasticity = 5e9, czEnergy = 50;
    int numberOfCores; // for job

    double indenterRadius=0.05;
    double indenterDepth=0;
    double indentationRate=0.001;
    double horizontalOffset=0;
    bool insertCZs = true;

    enum ConfinementType { bottomOnly=0, full=1, sides=2, frontAndBack=3 };

    ConfinementType confinement = ConfinementType::full;

    QString batchFileName;          // where this configuration is saved
    QString BatchName() const;      // same as batchFileName, without extension and path

    void Load(QString fileName);

    static QString ListToString(QList<double> &vec);
    static QString ToShortName(QString filePath);

    // GENERATION OF PY FILES
    void GeneratePythonScript();      // save .PY files into a fodler
};

#endif // BATCHCONFIGURATION_H
