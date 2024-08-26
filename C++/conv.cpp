#include <iostream>
#include <limits>  // For machine epsilon
#include <numeric> // For std::accumulate
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath> // For math functions like pow
#include <algorithm> // For std::max

using namespace std;

std::vector<double> readCSVArray(const std::string& filename) {
    std::vector<double> array1;
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Failed to open file." << std::endl;
        return array1;  // Return empty vector on failure
    }

    std::string line;
    if (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string value;

        // Read the comma-separated values into the array
        while (std::getline(ss, value, ',')) {
            array1.push_back(std::stod(value));
        }
    }

    file.close();
    return array1;
}
bool readCSV(const std::string& filename, std::vector<double>& x, std::vector<double>& y) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file." << std::endl;
        return false;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::stringstream lineStream(line);
        std::string col1, col2;

        if (std::getline(lineStream, col1, ',') && std::getline(lineStream, col2, ',')) {
            try {
                double x_val = std::stod(col1);  // Convert the first column to double
                double y_val = std::stod(col2);  // Convert the second column to double
                x.push_back(x_val);
                y.push_back(y_val);
            } catch (const std::invalid_argument& e) {
                std::cerr << "Invalid argument: " << e.what() << " for line: " << line << std::endl;
                return false;
            } catch (const std::out_of_range& e) {
                std::cerr << "Out of range: " << e.what() << " for line: " << line << std::endl;
                return false;
            }
        } else {
            std::cerr << "Error parsing line: " << line << std::endl;
            return false;
        }
    }

    file.close();
    return true;
}
void model(double Vp, std::vector<double> array1, double& RR, double& Ie, double& In, double& Ip, double& Io);
// Define a model function with an arbitrary number of parameters
int main() {

std::vector<double> array1 = readCSVArray("pa.csv");
// Array input

std::vector<double> V;
std::vector<double> I;
// Read the CSV file
if (!readCSV("in.csv", V, I)) {
std::cout << "Fail to find in.csv\n";
return 1; // Exit if reading failed
}
// Create and open a CSV file
std::ofstream csvFile("out.csv");
// Check if the file was opened successfully
if (csvFile.is_open())
{
    for(size_t i = 0; i < V.size(); ++i)
    {
     double rs,Ie,In,Ip,Io;
     model(V[i], array1,rs, Ie, In, Ip, Io);
     csvFile << V[i] << "," << Io <<"," << Ie <<"," << Ip <<"," << In <<"," << rs << "\n";
    }
    // Close the file
    csvFile.close();
    //std::cout << "Data written to data.csv successfully.\n";
    std::cout << " Te= " << array1[0] << " Tn= " << array1[1] << " Tp= " << array1[2] << " np= "  << array1[3]   << "\n";
    std::cout << " as= " << array1[4] << " Vs =" << array1[5] << " mode=" << array1[9] <<" rss=" << array1[6] << " Spap=" << array1[7] << " B=" << array1[8] << "\n";
    std::cout << " mr=" << array1[10] << "\n";

}
else
{
std::cerr << "Unable to open file.\n";
}

}
