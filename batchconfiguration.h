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


class BatchConfiguration
{
public:

    QString mshFileName;
    int dims = 2;

    double indenterRadius = 0.161925;
    double indentationDepth = 0.05;
    double horizontalOffset = 0.2;
    bool plasticity = true;

    double timeToRun = 0.2;
    int nFrames = 200;

    bool loadWithIndenter = true;   // if false -> static load

    constexpr static double blockHeight = 0.5;
    constexpr static double blockLength = 0.6;
    constexpr static double indentationRate = 0.2;

    constexpr static int numberOfCores = 12;
    constexpr static double YoungsModulus = 9e9;
    constexpr static double czsStrength = 4e5;
    constexpr static double czElasticity = 5e9;
    constexpr static double czEnergy = 50;

    bool insertCZs = true;


    QString batchFileName;          // where this configuration is saved
    QString BatchName() const;      // same as batchFileName, without extension and path

    void Load(QString fileName);

    static QString ListToString(QList<double> &vec);
    static QString ToShortName(QString filePath);

    // GENERATION OF PY FILES
    void GeneratePythonScript();      // save .PY files into a fodler

private:
    void GeneratePythonScript2D();
    void GeneratePythonScript3D();


};

#endif // BATCHCONFIGURATION_H
