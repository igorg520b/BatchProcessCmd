#include "batchconfiguration.h"

#include <QFileInfo>
#include <QTextStream>

#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <string>
#include <map>

#include <spdlog/spdlog.h>


QString BatchConfiguration::BatchName() const
{
    return ToShortName(batchFileName);
}


void BatchConfiguration::Load(QString fileName)
{
    QFile f(fileName);
    f.open(QIODevice::ReadOnly | QIODevice::Text);
    if(!f.exists()) qDebug() << "config file does not exist";
    QTextStream ts(&f);

    if (!f.exists())
    {
        qDebug() << "config file not open";
        throw std::runtime_error("");
    }

    QMap<QString, QString> kvp;
    QString line;
    while(!ts.atEnd())
    {
        line = ts.readLine();
        QStringList split = line.split("=");
        kvp.insert(split[0],split[1]);
    }
    f.close();

    if(kvp.contains("dim")) dim = kvp.value("dim").toInt();



    mshFileName = kvp.value("mesh");
    if(kvp.contains("czsStrength")) czsStrength = kvp.value("czStrength").toDouble();
    if(kvp.contains("YoungsModulus")) YoungsModulus = kvp.value("YoungsModulus").toDouble();

    if(kvp.contains("czElasticity")) czElasticity = kvp.value("czElasticity").toDouble();
    if(kvp.contains("czEnergy")) czEnergy = kvp.value("czEnergy").toDouble();

    numberOfCores = kvp.value("nCPUs").toInt();

    if(kvp.contains("indenterRadius")) indenterRadius = kvp.value("indenterRadius").toDouble();
    if(kvp.contains("indenterDepth")) indenterDepth = kvp.value("indenterDepth").toDouble();
    if(kvp.contains("indentationRate")) indentationRate = kvp.value("indentationRate").toDouble();
    if(kvp.contains("horizontalOffset")) horizontalOffset = kvp.value("horizontalOffset").toDouble();

    if(kvp.contains("insertCZs")) insertCZs = kvp.value("insertCZs")!="false";

    if(kvp.contains("confinement"))
    {
        QString s = kvp.value("confinement");
        if(s=="bottomOnly") confinement=ConfinementType::bottomOnly;
        else if(s=="frontAndBack") confinement=ConfinementType::frontAndBack;
        else if(s=="full") confinement=ConfinementType::full;
        else if(s=="sides") confinement=ConfinementType::sides;
    }
//{1: ^{0}d}|{2: ^{0}.3e}
    spdlog::info("test");
    spdlog::info("{0:<15} : {1}", "mshFileName", mshFileName.toStdString().c_str());
    spdlog::info("{0:<15} : {1}", "dim", dim);
    spdlog::info("{0:<15} : {1:03.2e}", "czsStrength", czsStrength);
    spdlog::info("{0:<15} : {1:03.2e}", "YoungsModulus", YoungsModulus);
    spdlog::info("{0:<15} : {1:03.2e}", "czElasticity", czElasticity);
    spdlog::info("{0:<15} : {1}", "czEnergy", czEnergy);
    spdlog::info("{0:<15} : {1}", "numberOfCores", numberOfCores);
    spdlog::info("{0:<15} : {1}", "indenterRadius", indenterRadius);
    spdlog::info("{0:<15} : {1}", "indenterDepth", indenterDepth);
    spdlog::info("{0:<15} : {1}", "indentationRate", indentationRate);
    spdlog::info("{0:<15} : {1}", "horizontalOffset", horizontalOffset);
    spdlog::info("{0:<15} : {1}", "insertCZs", insertCZs);
    spdlog::info("{0:<15} : {1}", "confinement", confinement);

    //
    batchFileName = fileName;
}



QString BatchConfiguration::ToShortName(QString filePath)
{
    QFileInfo fi(filePath);
    return fi.baseName();
}
/*
void BatchConfiguration::PrepareTable()
{
    qDebug() << "BatchConfiguration::PrepareTable";

}

*/


void BatchConfiguration::GeneratePythonScript()
{
    spdlog::info("BatchConfiguration::GeneratePythonScript()");

    // create directory for .py file
    QString tentativeDir = QDir::currentPath()+ "/" + BatchName();
    QDir dir(tentativeDir);
    if(!dir.exists()) dir.mkdir(tentativeDir);


    /*
    qDebug() << "BatchConfiguration::ProducePYFiles";
    icy::Mesh m;
    for(const TableEntry &te : qAsConst(tableEntries))
    {
        QString meshPath = "meshes/"+te.mshFileName;
        QFileInfo f(meshPath);
        if(!f.exists())
        {
            qDebug() << "mesh file not found" << meshPath;
            throw std::runtime_error("mesh file does not exist");
        }

        QString pyPath = QDir::currentPath()+ "/" + BatchName() + "/py/" + te.pyFileName;
        qDebug() << "loading mesh " << meshPath;
        m.LoadMSH(meshPath.toStdString(), !RHITA, insertCZs);
        m.RotateSample(te.rotationAngle);
        QString taskName = BatchName()+"_"+QString{"%1"}.arg(te.id,4,10,QLatin1Char('0'));
        m.ExportForAbaqus(pyPath.toStdString(), te.czStrength,taskName.toStdString(), BatchName().toStdString(),
                          YoungsModulus, czElasticity, czEnergy,
                          RHITA, indenterRadius, indenterDepth,indentationRate, horizontalOffset,numberOfCores,
                          (int)confinement);
    }
*/
}

