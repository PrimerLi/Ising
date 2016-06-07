#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
using namespace std;

double uniform()
{
    int randomNumber = rand();
    return float(randomNumber)/RAND_MAX;
}

double calMagnetization(const vector< vector<int> > &spin)
{
    double s = 0;
    int row = spin.size();
    int col = spin[0].size();
    for (int i = 0; i < row; ++i)
    {
        for (int j = 0; j < col; ++j)
        {
            s = s + spin[i][j];
        }
    }
    return s/float(row*col);
}

struct Neighbor
{
    int row, col;
    vector<int> left, right;
    vector<int> up, down;
    Neighbor(int row, int col)
    {
        this->row = row;
        this->col = col;
        left.push_back(row-1);
        for (int i = 1; i < row; ++i)
        {
            left.push_back(i-1);
        }
       for (int i = 0; i < row-1; ++i)
       {
           right.push_back(i+1);
       }
       right.push_back(0);
       up.push_back(col-1);
       for (int i = 1; i < col; ++i)
       {
            up.push_back(i-1);
       }
       for (int i = 0; i < col-1; ++i)
       {
           down.push_back(i+1);
       }
       down.push_back(0);
    }
    
    void print() const
    {
        cout << "Left: " << endl;
        for (int i = 0; i < left.size(); ++i)
        {
            cout << left[i] << "    ";
        }
        cout << endl << "Right: " << endl;
        for (int i = 0; i < right.size(); ++i)
        {
            cout << right[i] << "    ";
        }
        cout << endl << "Up: " << endl;
        for (int i = 0; i < up.size(); ++i)
        {
            cout << up[i] << "    ";
        }
        cout << endl << "Down: " << endl;
        for (int i = 0; i < down.size(); ++i)
        {
            cout << down[i] << "    ";
        }
        cout << endl;
    }
};

vector<int> generateNeighbors(const vector< vector<int> >&spin, int i, int j)
{
    vector<int> neighbors;
    int row = spin.size();
    int col = spin[0].size();
    if (i < 0 || i > row-1 || j < 0 || j > col - 1)
    {
        cout << "Wrong index. " << endl;
        exit(-1);
    }
    int left, right, up, down;
    Neighbor neighborTable(row, col);
    left = spin[neighborTable.left[i]][j];
    right = spin[neighborTable.right[i]][j];
    up = spin[i][neighborTable.up[j]];
    down = spin[i][neighborTable.down[j]];
    neighbors.push_back(left);
    neighbors.push_back(right);
    neighbors.push_back(up);
    neighbors.push_back(down);
    return neighbors;
}

double calEnergy(const vector< vector<int> > &spin, double h)
{
    double s = 0;
    int row = spin.size();
    int col = spin[0].size();
    Neighbor neighbors(row, col);
    for (int i = 0; i < row; ++i)
    {
        for (int j = 0; j < col; ++j)
        {
            int left, right;
            int up, down;
            left = spin[neighbors.left[i]][j];
            right = spin[neighbors.right[i]][j];
            up = spin[i][neighbors.up[j]];
            down = spin[i][neighbors.down[j]];
            s = s + spin[i][j]*(left + right + up + down);
        }
    }
    s = -0.5*s/float(row*col);
    double magnetization = calMagnetization(spin);
    s = s - h*magnetization;
    return s;
}

void createSpins(vector< vector<int> > &spin, int row, int col)
{
    for (int i = 0; i < row; ++i)
    {
        vector<int> temp;
        for (int j = 0; j < col; ++j)
        {
            double randomNumber = uniform();
            int sigma;
            if (randomNumber > 0.5)
                sigma = 1;
            else
                sigma = -1;
            temp.push_back(sigma);
        }
        spin.push_back(temp);
    }
}

void resetSpins(vector< vector<int> > &spin)
{
    int row = spin.size();
    int col = spin[0].size();
    for (int i = 0; i < row; ++i)
    {
        for (int j = 0; j < col; ++j)
        {
            int sigma;
            if (uniform() > 0.5)
                sigma = 1;
            else
                sigma = -1;
            spin[i][j] = sigma;
        }
    }
}

void printSpins(const vector< vector<int> > &spin)
{
    int row, col;
    row = spin.size();
    col = spin[0].size();
    for (int i = 0; i < row; ++i)
    {
        for (int j = 0; j < col; ++j)
        {
            cout << spin[i][j] << "    ";
        }
        cout << endl;
    }
    cout << endl;
}

int nint(double a)
{
    if (a <= 0)
        return -nint(-a);
    int int_part;
    int_part = int(a);
    double diff;
    diff = a-int_part;
    if (diff > 0.5)
        return int_part+1;
    else if (diff < 0.5)
        return int_part;
    else
    {
        double randomNumber = uniform();
        if (randomNumber > 0.5)
            return int_part;
        else
            return int_part + 1;
    }
}

void update(vector< vector<int> > &spin, double h, double beta)
{
    int row = spin.size();
    int col = spin[0].size();
    Neighbor neighbors(row, col);
    int randomRow, randomCol;
    randomRow = int(row*uniform());
    randomCol = int(col*uniform());
    bool accepted = false;
    int left, right, up, down;
    vector<int> myNeighbors = generateNeighbors(spin, randomRow, randomCol);
    left = myNeighbors[0];
    right = myNeighbors[1];
    up = myNeighbors[2];
    down = myNeighbors[3];
    double originalEnergy = -spin[randomRow][randomCol]*(left+right+up+down);
    double updatedEnergy = spin[randomRow][randomCol]*(left+right+up+down);
    originalEnergy = originalEnergy + (-h*spin[randomRow][randomCol]);
    updatedEnergy = updatedEnergy + (h*spin[randomRow][randomCol]);
    double energyDiff = updatedEnergy - originalEnergy;
    if (energyDiff > 0)
    {
        if (exp(-beta*energyDiff) > uniform())
        {
            accepted = true;
        }
    }
    else
    {
        accepted = true;
    }
    if (accepted)
    {
        spin[randomRow][randomCol] = -spin[randomRow][randomCol];
    }
}

void warmUp(vector< vector<int> > &spin, int warm, double h, double beta)
{
    for (int i = 0; i < warm; ++i)
    {
        update(spin, h, beta);
    }
}

double mean(const vector<double> &array)
{
    double s = 0;
    for (int i = 0; i < array.size(); ++i)
        s = s + array[i];
    return s/float(array.size());
}

double standardDeviation(const vector<double> &array)
{
    double mu = mean(array);
    double s = 0;
    int length = array.size();
    for (int i = 0; i < length; ++i)
        s = s + (array[i] - mu)*(array[i] - mu);
    s = s/float(length);
    return sqrt(s);
}

int main(int argc, char **argv)
{
    int row = 30;
    int col = 30;
    vector< vector<int> > spin;
    createSpins(spin, row, col);
    int warm = 8000;

    if (true)
    {
        double h1 = 3;
        double h2 = -3;
        int nh = 15;
        double deltah = (h2 - h1)/float(nh);
        double T;
        cout << "T = ";
        cin >> T;
        double beta = 1.0/T;
        vector<double> hvalues;
        int skip = 500;
        int numberOfMeasurements = 50;
        vector<double> energy, magnetization;
        for (int i = 0; i < nh; ++i)
        {
            double h = h1 + i*deltah;
            cout << h << endl;
            hvalues.push_back(h);
            warmUp(spin, warm, h, beta);
            vector<double> energyTemp, magnetizationTemp;
            for (int j = 0; j < numberOfMeasurements; ++j)
            {
                energyTemp.push_back(calEnergy(spin, h));
                magnetizationTemp.push_back(calMagnetization(spin));
                warmUp(spin, skip, h, beta);
            }
            energy.push_back(mean(energyTemp));
            magnetization.push_back(mean(magnetizationTemp));
        }
        for (int i = 0; i < nh; ++i)
        {
            double h = h2 - i*deltah;
            cout << h << endl;
            hvalues.push_back(h);
            warmUp(spin, warm, h, beta);
            vector<double> energyTemp, magnetizationTemp;
            for (int j = 0; j < numberOfMeasurements; ++j)
            {
                energyTemp.push_back(calEnergy(spin, h));
                magnetizationTemp.push_back(calMagnetization(spin));
                warmUp(spin, skip, h, beta);
            }
            energy.push_back(mean(energyTemp));
            magnetization.push_back(mean(magnetizationTemp));
        }
        ofstream ofile;
        string fileName = "energy_h_" + to_string(T) + ".txt";
        ofile.open(fileName);
        for (int i = 0; i < energy.size(); ++i)
        {
            ofile << hvalues[i] << "   " << energy[i] << endl;
        }
        ofile.close();
        fileName = "magnetization_h_" + to_string(T) + ".txt";
        ofile.open(fileName);
        for (int i = 0; i < magnetization.size(); ++i)
        {
            ofile << hvalues[i] << "    " << magnetization[i] << endl;
        }
        ofile.close();
    }
    
    if (false)
    {
        vector<double> energy;
        vector<double> magnetization;
        vector<double> specificHeat;
        vector<double> magneticSusceptibility;
        vector<double> temperature;
        double T1 = 6;
        double T2 = 0.001;
        double deltaT;
        int N = 20;
        int skip = 1000;
        int numberOfMeasurements = 100;
        double h = 0;
        deltaT = (T2 - T1)/float(N);
        for (int i = 0; i < N; ++i)
        {
            double T = T1 + i*deltaT;
            cout << T << endl;
            temperature.push_back(T);
            double beta = 1.0/T;
            //resetSpins(spin);
            warmUp(spin, warm, h, beta);
            vector<double> energyTemp, magnetizationTemp;
            for (int j = 0; j < numberOfMeasurements; ++j)
            {
                warmUp(spin, skip, h, beta);
                energyTemp.push_back(calEnergy(spin, h));
                magnetizationTemp.push_back(fabs(calMagnetization(spin)));
            }
            energy.push_back(mean(energyTemp));
            magnetization.push_back(mean(magnetizationTemp));
            specificHeat.push_back(standardDeviation(energyTemp));
            magneticSusceptibility.push_back(standardDeviation(magnetizationTemp));
        }
        ofstream ofile;
        ofile.open("energy.txt");
        for (int i = 0; i < energy.size(); ++i)
        {
            ofile << temperature[i] << "    " << energy[i] << endl;
        }
        ofile.close();
        ofile.open("magnetization.txt");
        for (int i = 0; i < magnetization.size(); ++i)
        {
            ofile << temperature[i] << "    " << magnetization[i] << endl;
        }
        ofile.close();
        ofile.open("specifiHeat.txt");
        for (int i = 0; i < specificHeat.size(); ++i)
        {
            ofile << temperature[i] << "    " << specificHeat[i] << endl;
        }
        ofile.close();
        ofile.open("magneticSusceptibility.txt");
        for (int i = 0; i < magneticSusceptibility.size(); ++i)
        {
            ofile << temperature[i] << "    " << magneticSusceptibility[i] << endl;
        }
        ofile.close();
    }
    return 0;
}
