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
#define RUN_FOREVER false  // set true if you want to run forever, but remember to set ITERATIONS_OF_GRASP to 0
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
vector<vector<int>> local_search(vector<vector<int>> solution, Transport transport, Customers customers, double **time_matrix, double *best);
vector<vector<int>> greedy_randomized(Transport transport, Customers customers, double **time_matrix, double *best);
double cost(vector<vector<int>> solution, Customers customers, double **time_matrix);

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

    if (!valid(customers, transport)) {
        // TODO edit save to file
        cout << "Error in data" << endl;
        return 1;
    }
    // create distance matrix
    // float **distanceMatrix;
    // distanceMatrix = new float *[customers.customers.size()];
    // for (int n = 0; n < customers.customers.size(); n++) {
    //     distanceMatrix[n] = new float[customers.customers.size()];
    // }
    // calculateDistanceMatrix(customers.customers, distanceMatrix);

    grasp(transport, customers);

    return 0;

    // customer home(0, 0, 0, 0, 0, 0, 0);
    // customers.push_back(home);
}
void chechIfDone(CodeExecutionCutoffTimer t) {
    if (t.isTimeUp()) {
        cout << "Time is up" << endl;
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
void saveToFile(vector<vector<int>> solution, double cost) {
    ofstream file("solution.txt");
    file << solution.size() << " " << setprecision(5) << fixed << cost << "\n";
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
            if (rand() % 2 == 0 && customers_copy.size() > 1)
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
    double best, local_best;
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
        cout << "iter: " << iter << endl;
        vector<vector<int>> instance = greedy_randomized(transport, customers, time_matrix, &best);
        vector<vector<int>> local_solution = local_search(instance, transport, customers, time_matrix, &local_best);
        // local_best = numeric_limits<double>::max();
        cout << "local best: " << local_best << endl;
        cout << "best: " << best << endl;
        if (local_best < best && local_best != -1) {
            best = local_best;
            saveToFile(local_solution, local_best);
        } else {
            saveToFile(instance, best);
        }
        chechIfDone(timer);
    }
    return best;
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

vector<vector<int>> local_search(vector<vector<int>> solution, Transport transport, Customers customers, double **time_matrix, double *best) {
    vector<vector<int>> new_solution;
    new_solution.assign(solution.begin(), solution.end());
    int iterations = 0;
    double best_cost = numeric_limits<double>::max();
    while (iterations++ < ITERATIONS_OF_LOCAL_SEARCH) {
        int i = rand() % new_solution.size();
        int j = rand() % new_solution.size();
        int k = rand() % new_solution[i].size();
        int l = rand() % new_solution[j].size();
        if (i == j && k == l) continue;
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
            *best = new_cost;
            best_cost = new_cost;
            iterations = 0;
        } else {
            new_solution[j][l] = new_solution[i][k];
            new_solution[i][k] = tmp;
        }
        chechIfDone(timer);
    }
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

    // double total_time = 0.0;
    // for (int veh = 0; veh < solution.size(); veh++) {
    //     double time = time_matrix[0][solution[veh][0]];
    //     if (time <= (double)customers.customers[solution[veh][0]].time_window_start) {
    //         time = (double)customers.customers[solution[veh][0]].time_window_start;
    //     }
    //     time += (double)customers.customers[solution[veh][0]].service_time;
    //     for (int pos = 1; pos < solution[veh].size(); pos++) {
    //         time += time_matrix[solution[veh][pos]][solution[veh][pos - 1]];
    //         if (time > (double)customers.customers[solution[veh][pos]].time_window_end) {
    //             return -1;
    //         }
    //         if (time <= (double)customers.customers[solution[veh][pos]].time_window_start) {
    //             time = (double)customers.customers[solution[veh][pos]].time_window_start;
    //         }
    //         time += (double)customers.customers[solution[veh][pos]].service_time;
    //     }
    //     time += time_matrix[solution[veh][solution[veh].size() - 1]][0];
    //     if (time > (double)customers.customers[0].time_window_end) {
    //         return -1;
    //     }
    //     total_time += time;
    // }
    // return total_time;
}

// int fitness(std::vector<std::vector<int>> individual, std::vector<customer> customers, int capacity) {
//     int total_fitness = 0;
//     // int temp_fitness = 0;
//     // vehicle dummy_vehicle(capacity);
//     // for (auto route : individual) {
//     //     for (int i = 1; i < customers.size(); i++) {
//     //         if (route[i] == 1) {
//     //             int t = time_to_arrive(dummy_vehicle, customers[i]);
//     //             if (t > customers[i].time_window_end) {
//     //                 return -1;
//     //             }
//     //             dummy_vehicle.capacity -= customers[i].demand;
//     //             dummy_vehicle.time += t + customers[i].service_time;
//     //             dummy_vehicle.x_cord = customers[i].x_cord;
//     //             dummy_vehicle.y_cord = customers[i].y_cord;
//     //             temp_fitness = dummy_vehicle.time;
//     //         }
//     //     }
//     //     dummy_vehicle.capacity = capacity;
//     //     dummy_vehicle.time = 0;
//     //     dummy_vehicle.x_cord = 0;
//     //     dummy_vehicle.y_cord = 0;
//     //     total_fitness += temp_fitness;
//     // }
//     return total_fitness;
// }

// bool is_acceptable(std::vector<std::vector<int>> individual, int number_of_vehicles, int number_of_customers) {
//     int check[number_of_customers] = {0};
//     for (int i = 0; i < number_of_vehicles; i++) {
//         for (int j = 0; j < number_of_customers; j++) {
//             if (check[i] + individual[i][j] > 1) {
//                 individual[i][j] = 0;
//             }
//             check[i] += individual[i][j];
//         }
//         if (check[i] != 1) {
//             return false;
//         }
//     }
//     return true;
// }

// int genetic_algorithm(std::vector<vehicle> vehicles, std::vector<customer> customers, int capacity, int number_of_individuals, float mutation_rate) {
//     int k = 0;
//     std::vector<std::vector<std::vector<int>>> population;
//     std::vector<int> fitnesses;
//     int min = -1;
//     for (int i = 0; i < number_of_individuals; i++) {
//         population.push_back(greedy_randomized(vehicles, customers));
//         fitnesses.push_back(fitness(population[i], customers, capacity));
//         min = fitnesses[i] < min || min == -1 ? fitnesses[i] : min;
//     }
//     while (k < 100) {
//         std::vector<std::vector<std::vector<int>>> new_population;
//         for (int i = 0; i < number_of_individuals; i++) {
//             fitnesses.clear();
//             // parent selection
//             int parent1 = rand() % number_of_individuals;
//             int parent2 = rand() % number_of_individuals;
//             while (parent1 == parent2) {
//                 parent2 = rand() % number_of_individuals;
//             }

//             int offspring1[vehicles.size()][customers.size() - 1], offspring2[vehicles.size()][customers.size() - 1];
//             std::vector<std::vector<int>> final_offspring;
//             final_offspring.resize(vehicles.size());

//             // crossover
//             for (int j = 0; j < vehicles.size(); j++) {
//                 final_offspring[j].resize(customers.size() - 1);
//                 for (int l = 0; l < customers.size() - 1; l++) {
//                     if (population[parent1][j][l] == population[parent2][j][l]) {
//                         offspring1[j][l] = population[parent1][j][l];
//                         offspring2[j][l] = population[parent1][j][l];
//                     } else {
//                         if (rand() % 2 == 0) {
//                             offspring1[j][l] = population[parent1][j][l];
//                             offspring2[j][l] = population[parent2][j][l];
//                         } else {
//                             offspring1[j][l] = population[parent2][j][l];
//                             offspring2[j][l] = population[parent1][j][l];
//                         }
//                     }
//                 }
//                 int decimal1 = 0, decimal2 = 0;
//                 for (int l = customers.size() - 2; l >= 0; l--) {
//                     decimal1 += pow(2, customers.size() - l - 2) * offspring1[j][l];
//                     decimal2 += pow(2, customers.size() - l - 2) * offspring2[j][l];
//                 }
//                 decimal1 = customers.size() - decimal1 - 1;
//                 decimal2 = customers.size() - decimal2 - 1;
//                 for (int l = customers.size() - 2; l >= 0; l--) {
//                     offspring1[j][l] = decimal1 % 2;
//                     decimal1 /= 2;
//                     offspring2[j][l] = decimal2 % 2;
//                     decimal2 /= 2;
//                     final_offspring[j][l] = (offspring1[j][l] & population[parent1][j][l]) | (offspring2[j][l] & population[parent2][j][l]);

//                     // mutation
//                     if (rand() % 100 < mutation_rate) {
//                         final_offspring[j][l] = ~final_offspring[j][l];
//                     }
//                 }
//                 int t = fitness(final_offspring, customers, capacity);
//                 if (is_acceptable(final_offspring, vehicles.size(), customers.size() - 1) && t != -1) {
//                     new_population[i].push_back(final_offspring[j]);
//                     fitnesses.push_back(t);
//                     min = t < min || min == -1 ? t : min;
//                 }
//             }
//         }
//         population.clear();
//         population = new_population;
//         k++;
//     }
//     return min;
// }
