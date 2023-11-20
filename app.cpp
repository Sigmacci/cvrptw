#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

typedef struct point2D {
    int x;
    int y;
} point2D;

typedef struct vehicle {
    int capacity;  // capacity of vehicle left
    int time;      // spent time
    point2D coordinates;

    vehicle(int capacity, point2D coordinates) {
        this->capacity = capacity;
        this->time = 0;
        this->coordinates = coordinates;
    };
} vehicle;
class Transport {
   public:
    int vehicle_cap;           // capacity of vehicle
    int time_limit;            // when dispatcher closes gates
    point2D coordinates;       // coordinates of dispatcher
    vector<vehicle> vehicles;  // storage of all vehicles

    Transport() {
        this->vehicle_cap = 0;
        this->time_limit = 0;
        this->coordinates = {0, 0};
    }
    /// @brief create new vehicle and add it to vehicles vector
    /// @return index of new vehicle
    int dispatchNewVehicle() {
        vehicle new_vehicle(this->vehicle_cap, this->coordinates);
        this->vehicles.push_back(new_vehicle);
        return this->vehicles.size() - 1;
    }
};

typedef struct customer {
    int id;
    point2D coordinates;
    int demand;
    int time_window_start;
    int time_window_end;
    int service_time;

    customer(int id, point2D coordinates, int demand, int time_window_start, int time_window_end, int service_time) {
        this->id = id;
        this->coordinates = coordinates;
        this->demand = demand;
        this->time_window_start = time_window_start;
        this->time_window_end = time_window_end;
        this->service_time = service_time;
    };
} customer;
class Customers {
   public:
    vector<customer> customers;
    Customers(){};
    /// @brief create new customer and add it to customers vector
    /// @return index of new customer
    int addCustomer(int id, point2D coordinates, int demand, int time_window_start, int time_window_end, int service_time) {
        customer new_customer(id, coordinates, demand, time_window_start, time_window_end, service_time);
        this->customers.push_back(new_customer);
        return this->customers.size() - 1;
    }
};

bool read_data_from_file(string path, Transport &transport, Customers &customers);
int time_to_arrive(vehicle v, customer r);

std::vector<std::vector<int>> greedy_randomized(std::vector<vehicle> vehicles, std::vector<customer> customers);
void heapify(std::vector<customer> customers, std::vector<double> coefficient, int n, int i);
void heapsort(std::vector<customer> customers, std::vector<double> coefficient, int n);

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
    // DEBUG: check tranport buffer
    // cout << transport.time_limit << endl
    //      << transport.vehicle_cap << endl
    //      << transport.vehicles.size() << endl
    //      << transport.x_cord_of_dispatcher << endl
    //      << transport.y_cord_of_dispatcher << endl;
    // DEBUG: check count of customers
    // cout << customers.customers.size() << endl;
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
    bool vehicle = false, customer = false;
    while (getline(file, line)) {
        if (customer) {
            // fix empty line after headers
            if (line.empty())
                continue;
            istringstream iss(line);
            int id, x_cord, y_cord, demand, time_window_start, time_window_end, service_time;
            iss >> id >> x_cord >> y_cord >> demand >> time_window_start >> time_window_end >> service_time;
            if (id == 0) {
                transport.coordinates = {x_cord, y_cord};
                transport.time_limit = time_window_end;
            } else {
                customers.addCustomer(id, {x_cord, y_cord}, demand, time_window_start, time_window_end, service_time);
            }
        } else if (vehicle) {
            istringstream iss(line);
            int tmp;
            iss >> tmp >> transport.vehicle_cap;
            vehicle = false;
        } else if (line.find("NUMBER") != string::npos) {
            vehicle = true;
        } else if (line.find("CUST") != string::npos) {
            customer = true;
        }
    }
    file.close();
    return !vehicle && customer;
}

/// @brief Calculate time from vehicle to traget customer
/// @return time neccessary to arrive to customer
float distance(vehicle v, customer r) {
    return sqrt(pow(v.coordinates.x - r.coordinates.x, 2) + pow(v.coordinates.y - r.coordinates.y, 2));
}

// void heapify(std::vector<customer> customers, std::vector<double> coefficient, int n, int i) {
//     int min = i;
//     int left = 2 * i + 1;
//     int right = 2 * i + 2;
//     if (left < n && coefficient[left] < coefficient[min])
//         min = left;
//     if (right < n && coefficient[right] < coefficient[min])
//         min = right;
//     if (min != i) {
//         std::swap(coefficient[i], coefficient[min]);
//         std::swap(customers[i], customers[min]);
//         heapify(customers, coefficient, n, min);
//     }
// }

// void heapsort(std::vector<customer> customers, std::vector<double> coefficient, int n) {
//     for (int i = n / 2 - 1; i >= 0; i--) {
//         heapify(customers, coefficient, n, i);
//     }
//     for (int i = 1; i < n; i++) {
//         std::swap(coefficient[0], coefficient[n - i]);
//         std::swap(customers[0], customers[n - i]);
//         heapify(customers, coefficient, n - i, 0);
//     }
// }

std::vector<std::vector<int>> greedy_randomized(Transport transport, Customers customers) {
    std::vector<std::vector<int>> result;
    std::vector<int> route = {0};
    route.resize(customers.size() - 1);
    std::vector<double> coefficients;
    std::vector<vehicle> vehicles_copy;
    std::vector<customer> customers_copy;
    std::copy(vehicles, vehicles_copy, vehicles.size());
    std::copy(customers, customers_copy, customers.size());

    int i = 1, j = 0, r = 0;

    while (i < customers_copy.size()) {
        for (int k = 0; k < customers_copy.size(); k++) {
            coefficients.push_back(customers_copy[k].demand / time_to_arrive(vehicles_copy[j], customers_copy[k]));
        }
        heapsort(customers_copy, coefficients, customers_copy.size());

        if (rand() % 2 == 0)
            r = (rand() % (customers_copy.size() - 1)) + 1;
        else
            r = i;

        int t = time_to_arrive(vehicles_copy[j], customers_copy[r]);
        if (vehicles_copy[j].capacity >= customers_copy[r].demand &&
            vehicles_copy[j].time + t <= customers_copy[r].time_window_end) {
            vehicles_copy[j].capacity -= customers_copy[r].demand;
            vehicles_copy[j].time += t + customers_copy[r].service_time;
            vehicles_copy[j].x_cord = customers_copy[r].x_cord;
            vehicles_copy[j].y_cord = customers_copy[r].y_cord;
            route[customers_copy[r].id - 1] = 1;
            customers_copy.erase(customers_copy.begin() + i);
            coefficients.clear();
        } else {
            if (i == customers_copy.size()) {
                i = 0;
                result.push_back(route);
                if (++j == vehicles_copy.size()) {
                    return {{}};
                }
                continue;
            }
        }
        if (r == i)
            i++;
    }
    return result;
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
