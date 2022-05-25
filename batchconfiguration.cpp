#include "batchconfiguration.h"

#include <QFileInfo>
#include <QTextStream>

#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>

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

    mshFileNames = ts.readLine().split(";");
    QStringList czsValues = ts.readLine().split(";");
    czsStrengths.clear();
    for(QString &str : czsValues)
        czsStrengths.push_back(str.toDouble());
    QStringList rotationAnglesTxt = ts.readLine().split(";");
    rotationAngles.clear();
    for(QString &str : rotationAnglesTxt)
        rotationAngles.push_back(str.toDouble());

    qDebug() << mshFileNames;
    qDebug() << czsStrengths;
    qDebug() << rotationAngles;

    batchFileName = fileName;
    qDebug() << "loaded";
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

                TableEntry te {count, str, s, a, pyFileName};
                tableEntries.push_back(te);
                std::cout << std::setw(3) << count;
                std::cout << "  | " << std::setw(15) << pyFileName.toStdString();
                std::cout << "  | " << std::setw(10) << str.toStdString();
                std::cout << "  | " << std::setw(6) << std::scientific << s;
                std::cout << "  | " << std::setw(6) << std::fixed << a << "  |\n";

                ofs << "\"" << pyFileName.toStdString() << "\"" << ",";
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
        m.LoadMSH(meshPath.toStdString());
        m.ExportForAbaqus(pyPath.toStdString(), te.czStrength);
    }

}

void BatchConfiguration::Convert_PY_to_INP(QString pathToAbaqus)
{
    qDebug() << "BatchConfiguration::Convert_PY_to_INP";

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