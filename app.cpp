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
#define RUN_FOREVER false               // set true if you want to run forever, but remember to set ITERATIONS_OF_GRASP to 0
#define ITERATIONS_OF_GRASP 280         // best tested value 280
#define RAND_PERCENTAGE 60              // best tested value 60 ,0-100 rand of greedy
#define ITERATIONS_OF_LOCAL_SEARCH 400  // best tested value 400
#define ITERATIONS_OF_POINT_SEARCH 500  // best tested value 500

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
        // return false;
        return elapsed_seconds.count() > targetTimeInSeconds;
    }
    double timeElapsed() {
        chrono::time_point<chrono::system_clock> end = chrono::system_clock::now();
        chrono::duration<double> elapsed_seconds = end - start;
        return elapsed_seconds.count();
    }
};

bool read_data_from_file(string path, Transport &transport, Customers &customers);
double time_to_arrive(int x1, int x2, int y1, int y2);
void saveToFile(vector<vector<int>> solution, double cost);
void chechIfDone(CodeExecutionCutoffTimer timer);
bool valid(Customers customers, Transport transport);
double grasp(Transport transport, Customers customers);
vector<vector<int>> local_search(vector<vector<int>> solution, Transport transport, Customers customers, double **time_matrix, double *best, double best_cost);
vector<vector<int>> greedy_randomized(Transport transport, Customers customers, double **time_matrix);
double cost(vector<vector<int>> solution, Customers customers, double **time_matrix);

CodeExecutionCutoffTimer timer(EXECUTION_TIME);

// int ttt = 3;  //+2
// void tmp(int n, double time) {
//     ifstream file("solution.txt", ios::in);
//     ofstream file2("zzzout.txt", ios::app);
//     file2 << n << " " << file.rdbuf() << " " << time << endl;
//     file.close();
//     file2.close();
// }

int main(int argc, char *argv[]) {
    // check if file name or rand % is given
    if (argc != 2) {
        cout << "app file" << endl;
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

    if (!valid(customers, transport)) {
        ofstream file("solution.txt");
        file << -1 << " " << setprecision(5) << fixed << 0.0 << "\n";
        file.close();
        cout << "Solution does not exist" << endl;
        return 1;
    }
    grasp(transport, customers);

    // for (int n = 27; n <= 1002; n += 25) {
    //     ttt = n;
    //     // define variables for storing data
    //     Transport transport;
    //     Customers customers;
    //     // read data from file
    //     if (!read_data_from_file(argv[1], transport, customers)) {
    //         cout << "Error in reading data from file" << endl;
    //         return 1;
    //     }

    //     if (!valid(customers, transport)) {
    //         ofstream file("solution.txt");
    //         file << -1 << " " << setprecision(5) << fixed << 0.0 << "\n";
    //         file.close();
    //         cout << "Solution does not exist" << endl;
    //         return 1;
    //     }
    //     cout << n << endl;
    //     chrono::time_point<chrono::system_clock> start = chrono::system_clock::now();
    //     grasp(transport, customers);
    //     chrono::time_point<chrono::system_clock> end = chrono::system_clock::now();
    //     chrono::duration<double> elapsed_seconds = end - start;
    //     tmp(n - 2, elapsed_seconds.count());
    // }
    cout << "Solution in file solution.txt" << endl;
    return 0;
}

void chechIfDone(CodeExecutionCutoffTimer t) {
    if (t.isTimeUp()) {
        cout << "Time is up" << endl
             << "Solution in file solution.txt" << endl;
        exit(0);
    }
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
                                   // while (getline(file, line) && ttt--) {  // customer data
        if (line.length() < 7)
            continue;
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
void saveToFile(vector<vector<int>> solution, double cost) {
    ofstream file("solution.txt");
    file << solution.size() << " " << setprecision(5) << fixed << cost << "\n";
    // file << solution.size() << " " << setprecision(5) << fixed << cost;
    for (int i = 0; i < solution.size(); i++) {
        for (int j = 0; j < solution[i].size(); j++) {
            file << solution[i][j] << " ";
        }
        file << "\n";
    }
    file.close();
}

bool valid(Customers customers, Transport transport) {
    for (customer customer : customers.customers) {
        if (customer.time_window_start > customer.time_window_end)
            return false;
        if (customer.time_window_start < 0)
            return false;
        // 0 -> dispatcher
        if (customer.time_window_start + customer.service_time >= customers.customers[0].time_window_end)
            return false;
        if (customer.service_time < 0)
            return false;
        if (customer.demand < 0)
            return false;
        if (customer.demand > transport.vehicle_cap)
            return false;
    }
    return true;
}

/// @brief Calculate time from vehicle to traget customer
/// @return time neccessary to arrive to customer
double time_to_arrive(int x1, int x2, int y1, int y2) {
    return sqrt(pow((double)x1 - (double)x2, 2) + pow((double)y1 - (double)y2, 2));
}

double time_to_arrive2(vehicle vehicle, customer customer) {
    double t = sqrt(
        pow(vehicle.x_cord - customer.x_cord, 2) +
        pow(vehicle.y_cord - customer.y_cord, 2));
    if (t < customer.time_window_start) {
        t = customer.time_window_start - vehicle.time;
    }
    return t;
}

vector<vector<int>> greedy_randomized(Transport transport, Customers customers, double **time_matrix, double *best) {
    vector<vector<int>> result;
    vector<int> route;
    vector<double> coefficients;
    vector<customer> customers_copy;
    customers_copy.assign(customers.customers.begin() + 1, customers.customers.end());

    int j = transport.dispatchNewVehicle(), r = 0;
    while (!customers_copy.empty()) {
        int prev_customer_id = 0;
        if (!route.empty()) {
            prev_customer_id = route.back();
        }
        double max = 0;
        int maxk = -1;
        for (int k = 0; k < customers_copy.size(); k++) {
            coefficients.push_back((double)customers_copy[k].demand / time_matrix[prev_customer_id][customers.customers[k].id]);
            if (max < coefficients.back()) {
                max = coefficients.back();
                maxk = k;
            };
        }

        int iter = 0;
        while (++iter < ITERATIONS_OF_POINT_SEARCH) {
            if (rand() % 101 <= RAND_PERCENTAGE && customers_copy.size() > 1)
                r = (rand() % (customers_copy.size() - 1)) + 1;
            else
                r = maxk;

            double tmpT = time_matrix[prev_customer_id][customers_copy[r].id] + transport.vehicles[j].time;
            if (tmpT <= (double)customers_copy[r].time_window_end && transport.vehicles[j].capacity >= customers_copy[r].demand && tmpT + (double)customers_copy[r].service_time + time_matrix[customers_copy[r].id][0] <= (double)customers.customers[0].time_window_end) {
                transport.vehicles[j].time = tmpT < (double)customers_copy[r].time_window_start ? (double)customers_copy[r].time_window_start : tmpT;
                transport.vehicles[j].time += (double)customers_copy[r].service_time;
                transport.vehicles[j].x_cord = customers_copy[r].x_cord;
                transport.vehicles[j].y_cord = customers_copy[r].y_cord;
                transport.vehicles[j].capacity -= customers_copy[r].demand;
                route.push_back(customers_copy[r].id);
                customers_copy.erase(customers_copy.begin() + r);
                coefficients.clear();
                iter = ITERATIONS_OF_POINT_SEARCH + 1;
            }
            chechIfDone(timer);
        }
        if (iter < ITERATIONS_OF_POINT_SEARCH + 1) {
            result.push_back(route);
            route.clear();
            j = transport.dispatchNewVehicle();
        }
    }
    result.push_back(route);
    *best = cost(result, customers, time_matrix);
    return result;
}

double grasp(Transport transport, Customers customers) {
    double root_cost, local_cost, best_cost = numeric_limits<double>::max();
    double **time_matrix = new double *[customers.customers.size()];
    for (int i = 0; i < customers.customers.size(); i++) {
        time_matrix[i] = new double[customers.customers.size()];
        for (int j = 0; j < i + 1; j++) {
            time_matrix[i][j] = time_to_arrive(customers.customers[i].x_cord, customers.customers[j].x_cord, customers.customers[i].y_cord, customers.customers[j].y_cord);
            time_matrix[j][i] = time_matrix[i][j];
        }
        chechIfDone(timer);
    }
    srand(time(NULL));
    int iter = 0;
    while (iter++ < ITERATIONS_OF_GRASP || RUN_FOREVER) {
        vector<vector<int>> instance = greedy_randomized(transport, customers, time_matrix, &root_cost);
        vector<vector<int>> local_solution = local_search(instance, transport, customers, time_matrix, &local_cost, root_cost);
        cout << "iteration: " << iter << " greedy_cost: " << setprecision(5) << fixed << root_cost << " local_cost: " << local_cost << endl;
        if (local_cost < best_cost && local_cost != -1) {
            best_cost = local_cost;
            saveToFile(local_solution, local_cost);
        } else if (root_cost < best_cost && root_cost != -1) {
            // if (root_cost < best_cost && root_cost != -1) {
            best_cost = root_cost;
            saveToFile(instance, root_cost);
        }
        chechIfDone(timer);
    }
    return best_cost;
}

bool isRouteValid(vector<int> route, int veh_cap, double **time_matrix, Customers customers) {
    double time = 0.0;
    int cap = 0;
    for (int n = 0; n < route.size(); n++) {
        time += time_matrix[route[n]][route[n - 1]];
        if (time > (double)customers.customers[route[n]].time_window_end) {
            return false;
        }
        if (time <= (double)customers.customers[route[n]].time_window_start) {
            time = (double)customers.customers[route[n]].time_window_start;
        }
        time += (double)customers.customers[route[n]].service_time;
        cap += customers.customers[route[n]].demand;
        if (cap > veh_cap) {
            return false;
        }
    }
    if (time + time_matrix[route[route.size() - 1]][0] > (double)customers.customers[0].time_window_end) {
        return false;
    }
    return true;
}

vector<vector<int>> local_search(vector<vector<int>> solution, Transport transport, Customers customers, double **time_matrix, double *best, double best_cost) {
    vector<vector<int>> new_solution;
    new_solution.assign(solution.begin(), solution.end());
    int iterations = 0;
    // double best_cost = numeric_limits<double>::max();
    while (iterations++ < ITERATIONS_OF_LOCAL_SEARCH) {
        int i = rand() % new_solution.size();
        int j = rand() % new_solution.size();
        int k = rand() % new_solution[i].size();
        int l = rand() % new_solution[j].size();
        if (i == j && k == l)
            continue;
        int tmp = new_solution[i][k];
        new_solution[i][k] = new_solution[j][l];
        new_solution[j][l] = tmp;

        if (!isRouteValid(new_solution[i], transport.vehicle_cap, time_matrix, customers) || !isRouteValid(new_solution[j], transport.vehicle_cap, time_matrix, customers)) {
            new_solution[j][l] = new_solution[i][k];
            new_solution[i][k] = tmp;
            continue;
        }

        double new_cost = cost(new_solution, customers, time_matrix);
        if (new_cost != -1 && new_cost < best_cost) {
            best_cost = new_cost;
            iterations = 0;
        } else {
            new_solution[j][l] = new_solution[i][k];
            new_solution[i][k] = tmp;
        }
        chechIfDone(timer);
    }
    *best = best_cost;
    // cout << best_cost << endl;
    return new_solution;
}
double distance(int x1, int y1, int x2, int y2) {
    return sqrt(pow((double)x1 - (double)x2, 2) + pow((double)y1 - (double)y2, 2));
}
double cost(vector<vector<int>> solution, Customers customers, double **time_matrix) {
    double total_time = 0.0;
    for (int veh = 0; veh < solution.size(); veh++) {
        double veh_time = 0.0;
        int prev_point = 0;
        for (int point = 0; point < solution[veh].size(); point++) {
            double t = time_matrix[prev_point][solution[veh][point]];
            // double t = distance(customers.customers[solution[veh][point]].x_cord, customers.customers[solution[veh][point]].y_cord, customers.customers[prev_point].x_cord, customers.customers[prev_point].y_cord);
            if (veh_time + t > (double)customers.customers[solution[veh][point]].time_window_end) {
                return -1;
            }
            if (veh_time + t <= (double)customers.customers[solution[veh][point]].time_window_start) {
                veh_time = (double)customers.customers[solution[veh][point]].time_window_start;
            } else {
                veh_time += t;
            }
            veh_time += (double)customers.customers[solution[veh][point]].service_time;
            prev_point = solution[veh][point];
        }
        veh_time += time_matrix[prev_point][0];
        // veh_time += distance(customers.customers[0].x_cord, customers.customers[0].y_cord, customers.customers[prev_point].x_cord, customers.customers[prev_point].y_cord);
        if (veh_time > (double)customers.customers[0].time_window_end) {
            return -1;
        }
        total_time += veh_time;
    }
    // cout << "total: " << total_time << endl;
    return total_time;
}
