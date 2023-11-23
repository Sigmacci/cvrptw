#include <time.h>

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#define EXECUTION_TIME 180.0

#define ITERATIONS_OF_LOCAL_SEARCH 100
#define ITERATIONS_OF_GRASP 1000
#define ITERATIONS_OF_POINT_SEARCH 100

using namespace std;

typedef struct vehicle {
    int capacity;  // capacity of vehicle left
    double time;   // spent time
    int x_cord;
    int y_cord;

    vehicle(int capacity, int x_cord, int y_cord) {
        this->capacity = capacity;
        this->time = 0;
        this->x_cord = x_cord;
        this->y_cord = y_cord;
    };
} vehicle;
class Transport {
   public:
    int vehicle_cap;           // capacity of vehicle
    int x_cord_of_dispatcher;  // x coordinate of dispatcher
    int y_cord_of_dispatcher;  // y coordinate of dispatcher
    vector<vehicle> vehicles;  // storage of all vehicles

    Transport() {
        this->vehicle_cap = 0;
        this->x_cord_of_dispatcher = 0;
        this->y_cord_of_dispatcher = 0;
    }
    /// @brief create new vehicle and add it to vehicles vector
    /// @return index of new vehicle
    int dispatchNewVehicle() {
        vehicle new_vehicle(this->vehicle_cap, this->x_cord_of_dispatcher, this->y_cord_of_dispatcher);
        this->vehicles.push_back(new_vehicle);
        return this->vehicles.size() - 1;
    }
};

typedef struct customer {
    int id;
    int x_cord;
    int y_cord;
    int demand;
    int time_window_start;
    int time_window_end;
    int service_time;

    customer(int id, int x_cord, int y_cord, int demand, int time_window_start, int time_window_end, int service_time) {
        this->id = id;
        this->x_cord = x_cord;
        this->y_cord = y_cord;
        this->demand = demand;
        this->time_window_start = time_window_start;
        this->time_window_end = time_window_end;
        this->service_time = service_time;
    };
} customer;
class Customers {
   public:
    vector<customer> customers;  // storage of all customers and point 0
    Customers(){};
    /// @brief create new customer and add it to customers vector
    /// @return index of new customer
    int addCustomer(int id, int x_cord, int y_cord, int demand, int time_window_start, int time_window_end, int service_time) {
        customer new_customer(id, x_cord, y_cord, demand, time_window_start, time_window_end, service_time);
        this->customers.push_back(new_customer);
        return this->customers.size() - 1;
    }
};

class CodeExecutionCutoffTimer {
   public:
    chrono::time_point<chrono::system_clock> start;
    double targetTimeInSeconds;
    CodeExecutionCutoffTimer(double targetTimeInSeconds) {
        this->targetTimeInSeconds = targetTimeInSeconds;
        this->start = chrono::system_clock::now();
    }
    bool isTimeUp() {
        chrono::time_point<chrono::system_clock> end = chrono::system_clock::now();
        chrono::duration<double> elapsed_seconds = end - start;
        return elapsed_seconds.count() > targetTimeInSeconds;
    }
};

bool read_data_from_file(string path, Transport &transport, Customers &customers);
double time_to_arrive(int x1, int x2, int y1, int y2);
void saveToFile(vector<vector<int>> solution, double cost);
void chechIfDone(CodeExecutionCutoffTimer timer);
bool valid(Customers customers, Transport transport);
double grasp(Transport transport, Customers customers);
double local_search(vector<vector<int>> solution, Transport transport, Customers customers, double **time_matrix);
vector<vector<int>> greedy_randomized(Transport transport, Customers customers, double **time_matrix);
double cost(vector<vector<int>> solution, Customers customers, double **time_matrix);
void readSolution(Customers customers);
CodeExecutionCutoffTimer timer(EXECUTION_TIME);

int main(int argc, char *argv[]) {
    // check if file name is given
    if (argc != 2) {
        cout << "Please enter file name" << endl;
        return 1;
    }
    // define variables for storing data
    Transport transport;
    Customers customers;
    // read data from file
    if (!read_data_from_file(argv[1], transport, customers)) {
        cout << "Error in reading data from file" << endl;
        return 1;
    }
    readSolution(customers);

    return 0;

    // customer home(0, 0, 0, 0, 0, 0, 0);
    // customers.push_back(home);
}

/// @brief Read neccecery data from file and store them in vectors
/// @param path name and path of file
/// @return true->success false->fail
bool read_data_from_file(string path, Transport &transport, Customers &customers) {
    // open file
    ifstream file(path);
    if (!file.is_open()) {
        return false;
    }
    string line;
    getline(file, line);  // name
    getline(file, line);  // empty
    getline(file, line);  // VEHICLE
    getline(file, line);  // HEADERS
    getline(file, line);  // vehicle data
    istringstream iss(line);
    int tmp;
    iss >> tmp >> transport.vehicle_cap;
    getline(file, line);           // empty
    getline(file, line);           // CUSTOMER
    getline(file, line);           // HEADERS
    while (getline(file, line)) {  // customer data
        if (line.length() < 7) continue;
        istringstream iss(line);
        int id, x_cord, y_cord, demand, time_window_start, time_window_end, service_time;
        iss >> id >> x_cord >> y_cord >> demand >> time_window_start >> time_window_end >> service_time;
        if (id == 0) {
            transport.x_cord_of_dispatcher = x_cord;
            transport.y_cord_of_dispatcher = y_cord;
        }
        customers.addCustomer(id, x_cord, y_cord, demand, time_window_start, time_window_end, service_time);
    }
    file.close();
    return true;
}
double distance(int x1, int y1, int x2, int y2) {
    return sqrt(pow((double)x1 - (double)x2, 2) + pow((double)y1 - (double)y2, 2));
}
void readSolution(Customers customers) {
    ifstream file("solution.txt");
    string line;
    getline(file, line);  // name
    istringstream iss(line);
    int tmp, prev = 0;
    double cost, calc_cost = 0.0;
    iss >> tmp >> cost;
    int linen = 0;
    while (getline(file, line)) {
        linen++;
        double time = 0.0;
        istringstream iss(line);
        while (iss >> tmp) {
            double t = distance(customers.customers[tmp].x_cord, customers.customers[tmp].y_cord, customers.customers[prev].x_cord, customers.customers[prev].y_cord);
            if (time + t <= (double)customers.customers[tmp].time_window_start) {
                time = (double)customers.customers[tmp].time_window_start;
            } else {
                time += t;
            }
            time += (double)customers.customers[tmp].service_time;
            prev = tmp;
        }
        time += distance(customers.customers[prev].x_cord, customers.customers[prev].y_cord, customers.customers[0].x_cord, customers.customers[0].y_cord);
        calc_cost += time;
        prev = 0;
    }
    cout << setprecision(6) << fixed << "cost: " << cost << " calc_cost: " << calc_cost << " line:" << linen << endl;
}
