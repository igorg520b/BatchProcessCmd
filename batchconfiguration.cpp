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

BatchConfiguration::BatchConfiguration()
{
    rotationAngles.push_back(0);
    czsStrengths.push_back(1.e5);
    czsStrengths.push_back(2.e5);
    czsStrengths.push_back(3.e5);
    mshFileNames.push_back("test1");
    mshFileNames.push_back("test2");
}

QString BatchConfiguration::BatchName() const
{
    return ToShortName(batchFileName);
}


void BatchConfiguration::Load(QString fileName)
{
    qDebug() << "BatchConfiguration::Load";

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

    czsStrengths.clear();
    rotationAngles.clear();

    mshFileNames = kvp.value("meshes").split(";");

    for(QString &str : kvp.value("czStrengths").split(";"))
        czsStrengths.push_back(str.toDouble());

    for(QString &str : kvp.value("rotationAngles").split(";"))
        rotationAngles.push_back(str.toDouble());

    pathToAbaqus = kvp.value("AbaqusPath");
    numberOfCores = kvp.value("nCPUs").toInt();
    YoungsModulus = kvp.value("YoungsModulus").toDouble();

    czElasticity = kvp.value("czElasticity").toDouble();
    czEnergy = kvp.value("czEnergy").toDouble();

    if(kvp.contains("RHITA")) RHITA = kvp.value("RHITA")=="y";
    if(kvp.contains("indenterRadius")) indenterRadius = kvp.value("indenterRadius").toDouble();
    if(kvp.contains("indenterDepth")) indenterDepth = kvp.value("indenterDepth").toDouble();
    if(kvp.contains("indentationRate")) indentationRate = kvp.value("indentationRate").toDouble();
    if(kvp.contains("horizontalOffset")) horizontalOffset = kvp.value("horizontalOffset").toDouble();

    qDebug() << "mshFileNames" << mshFileNames;
    qDebug() << "czsStrengths" << czsStrengths;
    qDebug() << "rotationAngles" << rotationAngles;
    qDebug() << "pathToAbaqus" << pathToAbaqus;
    qDebug() << "numberOfCores" << numberOfCores;
    qDebug() << "YoungsModulus" << YoungsModulus;
    qDebug() << "czElasticity" << czElasticity;
    qDebug() << "czEnergy" << czEnergy;
    qDebug() << "RHITA" << RHITA;
    qDebug() << "indenterRadius" << indenterRadius;
    qDebug() << "indenterDepth" << indenterDepth;
    qDebug() << "indentationRate" << indentationRate;
    qDebug() << "horizontalOffset" << horizontalOffset;

    //
    batchFileName = fileName;
}

void BatchConfiguration::ParseValueList(QString CSV, QList<double> &vect)
{
    vect.clear();
    std::stringstream ss(CSV.toStdString());

    for (double d; ss >> d;) {
        vect.push_back(d);
        if (ss.peek() == ',')
            ss.ignore();
    }
}

QString BatchConfiguration::ToShortName(QString filePath)
{
    QFileInfo fi(filePath);
    return fi.baseName();
}

void BatchConfiguration::PrepareTable()
{
    qDebug() << "BatchConfiguration::PrepareTable";

    QString tentativeDir = QDir::currentPath()+ "/" + BatchName();
    QDir dir(tentativeDir);
    if(!dir.exists())
    {
        dir.mkdir(tentativeDir);
        dir.mkdir(tentativeDir+"/py");
        dir.mkdir(tentativeDir+"/inp");
        dir.mkdir(tentativeDir+"/scratch");
    }
    QString outputCSV = dir.path() + "/" + BatchName() + ".csv";
    std::ofstream ofs(outputCSV.toStdString(), std::ios_base::trunc|std::ios_base::out);


    int count = 1;
    tableEntries.clear();
    for(QString &str : mshFileNames)
        for(const double &a : qAsConst(rotationAngles))
            for(const double &s : qAsConst(czsStrengths))
            {
                QString pyFileName = BatchName()+"_"+QString{"%1"}.arg(count,4,10,QLatin1Char('0'))+".py";
                QString taskName = BatchName()+"_"+QString{"%1"}.arg(count,4,10,QLatin1Char('0'));

                TableEntry te {count, str, s, a, pyFileName};
                tableEntries.push_back(te);
                std::cout << std::setw(3) << count;
                std::cout << "  | " << std::setw(15) << pyFileName.toStdString();
                std::cout << "  | " << std::setw(10) << str.toStdString();
                std::cout << "  | " << std::setw(6) << std::scientific << s;
                std::cout << "  | " << std::setw(6) << std::fixed << a << "  |\n";

                ofs << "\"" << taskName.toStdString() << "\"" << ",";
                ofs << "\"" << str.toStdString() << "\"" <<",";
                ofs << std::setw(6) << std::scientific << s << ",";
                ofs <<  std::setw(6) << std::fixed << a << "\n";
                count++;
            }
    ofs.close();
}

void BatchConfiguration::ProducePYFiles()
{
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
        m.LoadMSH(meshPath.toStdString(), !RHITA);
        m.RotateSample(te.rotationAngle);
        QString taskName = BatchName()+"_"+QString{"%1"}.arg(te.id,4,10,QLatin1Char('0'));
        m.ExportForAbaqus(pyPath.toStdString(), te.czStrength,taskName.toStdString(), BatchName().toStdString(),
                          YoungsModulus, czElasticity, czEnergy,
                          RHITA, indenterRadius, indenterDepth,indentationRate, horizontalOffset,numberOfCores);
    }

}

void BatchConfiguration::Convert_PY_to_INP()
{
    qDebug() << "BatchConfiguration::Convert_PY_to_INP";
    QString inpPath = "chdir " + QDir::currentPath()+ "\\" + BatchName() + "\\inp";
    std::system(inpPath.toStdString().c_str());

    QString batFileName = QDir::currentPath()+ "\\" + BatchName() + "\\inp\\" + BatchName() + ".bat";

    std::ofstream ofs(batFileName.toStdString(), std::ios_base::trunc|std::ios_base::out);

    for(const TableEntry &te : qAsConst(tableEntries))
    {
        QString pyPath = QDir::currentPath()+ "/" + BatchName() + "/py/" + te.pyFileName;
        QString command = pathToAbaqus + " cae noGUI=\"" + pyPath + "\"";
        std::system(command.toStdString().c_str());

        ofs << pathToAbaqus.toStdString() << " job=";
        QString taskName = BatchName()+"_"+QString{"%1"}.arg(te.id,4,10,QLatin1Char('0'));
        ofs << taskName.toStdString();
        ofs << " cpus=" << this->numberOfCores;
        ofs << " interactive\n";
    }

    ofs.close();
}


QString BatchConfiguration::ListToString(QList<double> &vec)
{
    QString r;
    for(int i=0;i<vec.size();i++)
    {
        r.append(QString{"%1"}.arg(vec[i],0,'e'));
        if(i!=vec.size()-1) r.append(',');
    }
    return r;
}
