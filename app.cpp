#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#define ITERATIONS_OF_LOCAL_SEARCH 100
#define ITERATIONS_OF_GRASP 100

using namespace std;

typedef struct vehicle {
    int capacity;  // capacity of vehicle left
    int time;      // spent time
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

bool read_data_from_file(string path, Transport &transport, Customers &customers);
float time_to_arrive(int x1, int x2, int y1, int y2);
void saveToFile(vector<vector<int>> solution);

bool valid(Customers customers, Transport transport);
int local_search(vector<vector<int>> solution, Transport transport, Customers customers, int **time_matrix);
int cost(vector<vector<int>> solution, Customers customers, int **time_matrix);
int grasp(Transport transport, Customers customers);

vector<vector<int>> greedy_randomized(Transport transport, Customers customers, int **time_matrix);
void heapify(vector<customer> customers, vector<double> coefficient, int n, int i);
void heapsort(vector<customer> customers, vector<double> coefficient, int n);

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

bool valid(Customers customers, Transport transport) {
    for (auto customer : customers.customers) {
        if (customer.time_window_start > customer.time_window_end)
            return false;
        if (customer.time_window_start < 0)
            return false;
        // 0 -> dispatcher
        if (customer.time_window_start + customer.service_time > customers.customers[0].time_window_end)
            return false;
        if (customer.time_window_end > customers.customers[0].time_window_end)
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
float time_to_arrive(int x1, int x2, int y1, int y2) {
    return sqrt(
        pow(x1 - x2, 2) +
        pow(y1 - y2, 2));
}

float time_to_arrive2(vehicle vehicle, customer customer) {
    float t = sqrt(
        pow(vehicle.x_cord - customer.x_cord, 2) +
        pow(vehicle.y_cord - customer.y_cord, 2));
    if (t < customer.time_window_start) {
        t = customer.time_window_start - vehicle.time;
    }
    return t;
}

/// @brief heapify function for heapsort
/// @param customers
/// @param coefficient
/// @param n
/// @param i
void heapify(vector<customer> customers, vector<double> coefficient, int n, int i) {
    int min = i;
    int left = 2 * i + 1;
    int right = 2 * i + 2;
    if (left < n && coefficient[left] < coefficient[min])
        min = left;
    if (right < n && coefficient[right] < coefficient[min])
        min = right;
    if (min != i) {
        swap(coefficient[i], coefficient[min]);
        swap(customers[i], customers[min]);
        heapify(customers, coefficient, n, min);
    }
}

/// @brief heapsort function
/// @param customers
/// @param coefficient
/// @param n
void heapsort(vector<customer> customers, vector<double> coefficient, int n) {
    for (int i = n / 2 - 1; i >= 0; i--) {
        heapify(customers, coefficient, n, i);
    }
    for (int i = 1; i < n; i++) {
        swap(coefficient[0], coefficient[n - i]);
        swap(customers[0], customers[n - i]);
        heapify(customers, coefficient, n - i, 0);
    }
}

vector<vector<int>> greedy_randomized(Transport transport, Customers customers, int **time_matrix) {
    vector<vector<int>> result;
    vector<int> route;
    // route.resize(customers.customers.size() - 1);
    vector<double> coefficients;
    vector<customer> customers_copy;
    customers_copy.assign(customers.customers.begin() + 1, customers.customers.end());

    int i = 1, j = transport.dispatchNewVehicle(), r = 0;
    while (!customers_copy.empty()) {
        int prev_customer_id = 0;
        if (!route.empty()) {
            prev_customer_id = route.back();
        }
        double max = 0;
        int maxk = -1;
        for (int k = 0; k < customers_copy.size(); k++) {
            coefficients.push_back((float)customers_copy[k].demand / (float)time_matrix[prev_customer_id][customers.customers[k].id + 1]);
            if (max < coefficients.back()) {
                max = coefficients.back();
                maxk = k;
            };
        }
        // heapsort(customers_copy, coefficients, customers_copy.size());

        if (rand() % 2 == 0 && customers_copy.size() > 1)
            r = (rand() % (customers_copy.size() - 1)) + 1;
        else
            r = maxk;

        int t = time_matrix[prev_customer_id][customers_copy[r].id];
        if (transport.vehicles[j].time + t < customers_copy[r].time_window_start) {
            t = customers_copy[r].time_window_start - transport.vehicles[j].time;
        }
        int ret_to_dep = time_matrix[customers_copy[r].id][0];
        if (transport.vehicles[j].capacity >= customers_copy[r].demand &&
            transport.vehicles[j].time + t <= customers_copy[r].time_window_end &&
            transport.vehicles[j].time + t + customers_copy[r].service_time + ret_to_dep <= customers.customers[0].time_window_end) {
            transport.vehicles[j].capacity -= customers_copy[r].demand;
            transport.vehicles[j].time += t + customers_copy[r].service_time;
            transport.vehicles[j].x_cord = customers_copy[r].x_cord;
            transport.vehicles[j].y_cord = customers_copy[r].y_cord;
            route.push_back(customers_copy[r].id);
            customers_copy.erase(customers_copy.begin() + r);
            coefficients.clear();
        } else {
            i = 0;
            result.push_back(route);
            route.clear();
            j = transport.dispatchNewVehicle();
            continue;
        }
        if (r == i) {
            i++;
        }
    }
    result.push_back(route);
    return result;
}

int grasp(Transport transport, Customers customers) {
    int best = numeric_limits<int>::max();
    int **time_matrix = new int *[customers.customers.size()];
    for (int i = 0; i < customers.customers.size(); i++) {
        time_matrix[i] = new int[customers.customers.size()];
        for (int j = 0; j < i + 1; j++) {
            time_matrix[i][j] = time_to_arrive(customers.customers[i].x_cord, customers.customers[j].x_cord, customers.customers[i].y_cord, customers.customers[j].y_cord);
            time_matrix[j][i] = time_matrix[i][j];
        }
    }
    int iter = 0, local_best = numeric_limits<int>::max();
    while (iter++ < ITERATIONS_OF_GRASP) {
        vector<vector<int>> instance = greedy_randomized(transport, customers, time_matrix);
        local_best = local_search(instance, transport, customers, time_matrix);
        if (local_best < best) {
            best = local_best;
            saveToFile(instance);
        }
    }
    cout << "UMM" << best << endl;

    return best;
}

int local_search(vector<vector<int>> solution, Transport transport, Customers customers, int **time_matrix) {
    vector<vector<int>> new_solution;
    new_solution.assign(solution.begin(), solution.end());
    int best_cost = cost(solution, customers, time_matrix);
    int iterations = 0;
    while (iterations++ < ITERATIONS_OF_LOCAL_SEARCH) {
        int i = rand() % new_solution.size();
        int j = rand() % new_solution.size();
        // while (i == j) {
        //     j = rand() % new_solution.size();
        // }

        int k = rand() % new_solution[i].size();
        int l = rand() % new_solution[j].size();
        swap(new_solution[i][k], new_solution[j][l]);
        int new_cost = cost(new_solution, customers, time_matrix);
        if (new_cost != -1 && new_cost < best_cost) {
            best_cost = new_cost;
            solution.clear();
            solution.assign(new_solution.begin(), new_solution.end());
            iterations = 0;
        } else {
            swap(new_solution[i][k], new_solution[j][l]);
        }
    }

    return best_cost;
}

int cost(vector<vector<int>> solution, Customers customers, int **time_matrix) {
    int total_time = 0;
    for (int veh = 0; veh < solution.size(); veh++) {
        int time = time_matrix[0][solution[veh][0]];
        if (time < customers.customers[solution[veh][0]].time_window_start) {
            time = customers.customers[solution[veh][0]].time_window_start;
        }
        time += customers.customers[solution[veh][0]].service_time;
        for (int pos = 1; pos < solution[veh].size(); pos++) {
            time += time_matrix[solution[veh][pos]][solution[veh][pos - 1]];
            if (time > customers.customers[solution[veh][pos]].time_window_end) {
                return -1;
            }
            if (time < customers.customers[solution[veh][pos]].time_window_start) {
                time = customers.customers[solution[veh][pos]].time_window_start;
            }
            time += customers.customers[solution[veh][pos]].service_time;
        }
        time += time_matrix[solution[veh][solution[veh].size() - 1]][0];
        if (time > customers.customers[0].time_window_end) {
            return -1;
        }
        total_time += time;
    }
    return total_time;
    //     int total_cost = 0;
    // for (int i = 0; i < solution.size(); i++) {
    //     vehicle dummy_vehicle(transport.vehicle_cap, transport.x_cord_of_dispatcher, transport.y_cord_of_dispatcher);
    //     for (int j = 0; j < solution[i].size(); j++) {
    //         int t = time_to_arrive2(dummy_vehicle, customers.customers[solution[i][j]]);
    //         if (t > customers.customers[solution[i][j]].time_window_end) {
    //             return -1;
    //         }
    //         dummy_vehicle.capacity -= customers.customers[i].demand;
    //         dummy_vehicle.time += t + customers.customers[i].service_time;
    //         dummy_vehicle.x_cord = customers.customers[i].x_cord;
    //         dummy_vehicle.y_cord = customers.customers[i].y_cord;
    //     }
    //     int t = time_to_arrive2(dummy_vehicle, customers.customers[0]);
    //     if (t > customers.customers[0].time_window_end) {
    //         return -1;
    //     }
    //     total_cost += dummy_vehicle.time + t;
    // }
    // return total_cost;
}

void saveToFile(vector<vector<int>> solution) {
    ofstream file("solution.txt");
    for (int i = 0; i < solution.size(); i++) {
        file << "Route " << i << ": ";
        for (int j = 0; j < solution[i].size(); j++) {
            file << solution[i][j] << " ";
        }
        file << endl;
    }
    file.close();
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
