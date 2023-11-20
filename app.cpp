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
class Transport
{
public:
    int vehicle_cap;          // capacity of vehicle
    int time;                 // when dispatcher closes gates
    int x_cord_of_dispatcher; // x coordinate of dispatcher
    int y_cord_of_dispatcher; // y coordinate of dispatcher
    vector<vehicle> vehicles; // storage of all vehicles

    Transport() {
        this->vehicle_cap = 0;
        this->time_limit = 0;
        this->coordinates = {0, 0};
    }
    Transport(){}; // temporary constructor

    Transport(int vehicle_cap, int x_cord_of_dispatcher, int y_cord_of_dispatcher)
    {
        this->vehicle_cap = vehicle_cap;
        this->time = 0;
        this->x_cord_of_dispatcher = x_cord_of_dispatcher;
        this->y_cord_of_dispatcher = y_cord_of_dispatcher;
    };
    /// @brief create new vehicle and add it to vehicles vector
    /// @return index of new vehicle
    int dispatchNewVehicle() {
        vehicle new_vehicle(this->vehicle_cap, this->coordinates);
    int dispatchNewVehicle()
    {
        vehicle new_vehicle(this->vehicle_cap, this->x_cord_of_dispatcher, this->y_cord_of_dispatcher);
        this->vehicles.push_back(new_vehicle);
        return this->vehicles.size() - 1;
    }
};

typedef struct customer
{
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
    vector<customer> customers;  // storage of all customers and point 0
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
float distance(point2D a, point2D b);
void calculateDistanceMatrix(vector<customer> points, float **distanceMatrix);

vector<vector<int>> greedy_randomized(Transport transport, Customers customers);
void heapify(vector<customer> customers, vector<double> coefficient, int n, int i);
void heapsort(vector<customer> customers, vector<double> coefficient, int n);

int main(int argc, char *argv[])
{
    // check if file name is given
    if (argc != 2)
    {
        cout << "Please enter file name" << endl;
        return 1;
    }
    // define variables for storing data
    Transport transport;
    Customers customers;
    // read data from file
    if (!read_data_from_file(argv[1], transport, customers))
    {
        cout << "Error in reading data from file" << endl;
        return 1;
    }
    // create distance matrix
    float **distanceMatrix;
    distanceMatrix = new float *[customers.customers.size()];
    for (int n = 0; n < customers.customers.size(); n++) {
        distanceMatrix[n] = new float[customers.customers.size()];
    }
    calculateDistanceMatrix(customers.customers, distanceMatrix);

    return 0;

    // customer home(0, 0, 0, 0, 0, 0, 0);
    // customers.push_back(home);
}

/// @brief Read neccecery data from file and store them in vectors
/// @param path name and path of file
/// @return true->success false->fail
bool read_data_from_file(string path, Transport &transport, Customers &customers, vector<customer> &points) {
    // open file
    ifstream file(path);
    if (!file.is_open())
    {
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
    getline(file, line);            // empty
    getline(file, line);            // CUSTOMER
    getline(file, line);            // HEADERS
    while (getline(file, line)) {   // customer data
        if (line == " ") continue;  // fix empty line after headers
        if (line.empty()) break;    // end of file
        istringstream iss(line);
        int id, x_cord, y_cord, demand, time_window_start, time_window_end, service_time;
        iss >> id >> x_cord >> y_cord >> demand >> time_window_start >> time_window_end >> service_time;
        if (id == 0) {
            transport.coordinates = {x_cord, y_cord};
            transport.time_limit = time_window_end;
        }
        customers.addCustomer(id, {x_cord, y_cord}, demand, time_window_start, time_window_end, service_time);
    }
    file.close();
    return true;
}

bool valid(Customers customers, Transport transport)
{
    for (auto customer : customers.customers)
    {
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

bool valid(Customers customers, Transport transport)
{
    for (auto customer : customers.customers)
    {
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
int time_to_arrive(vehicle v, customer r)
{
    float t = sqrt(
        pow(v.x_cord - r.x_cord, 2) +
        pow(v.y_cord - r.y_cord, 2));

    if (v.time + t < r.time_window_start)
        return r.time_window_start - v.time;

    return t;
}

/// @brief heapify function for heapsort
/// @param customers 
/// @param coefficient 
/// @param n 
/// @param i 
void heapify(vector<customer> customers, vector<double> coefficient, int n, int i)
{
    int min = i;
    int left = 2 * i + 1;
    int right = 2 * i + 2;
    if (left < n && coefficient[left] < coefficient[min])
        min = left;
    if (right < n && coefficient[right] < coefficient[min])
        min = right;
    if (min != i)
    {
        swap(coefficient[i], coefficient[min]);
        swap(customers[i], customers[min]);
        heapify(customers, coefficient, n, min);
    }
}

/// @brief heapsort function
/// @param customers 
/// @param coefficient 
/// @param n 
void heapsort(vector<customer> customers, vector<double> coefficient, int n)
{
    for (int i = n / 2 - 1; i >= 0; i--)
    {
        heapify(customers, coefficient, n, i);
    }
    for (int i = 1; i < n; i++)
    {
        swap(coefficient[0], coefficient[n - i]);
        swap(customers[0], customers[n - i]);
        heapify(customers, coefficient, n - i, 0);
    }
}

vector<vector<int>> greedy_randomized(Transport transport, Customers customers)
{
    vector<vector<int>> result;
    vector<int> route;
    // route.resize(customers.customers.size() - 1);
    vector<double> coefficients;
    vector<vehicle> vehicles_copy;
    vector<customer> customers_copy;
    copy(transport.vehicles, vehicles_copy, transport.vehicles.size());
    copy(customers.customers, customers_copy, customers.customers.size());

//     int i = 1, j = 0, r = 0;

    while (i < customers_copy.size())
    {
        for (int k = 0; k < customers_copy.size(); k++)
        {
            coefficients.push_back(customers_copy[k].demand / time_to_arrive(vehicles_copy[j], customers_copy[k]));
        }
        heapsort(customers_copy, coefficients, customers_copy.size());

//         if (rand() % 2 == 0)
//             r = (rand() % (customers_copy.size() - 1)) + 1;
//         else
//             r = i;

        int t = time_to_arrive(vehicles_copy[j], customers_copy[r]);
        int ret_to_dep = time_to_arrive(vehicles_copy[j], customers_copy[0]);
        if (vehicles_copy[j].capacity >= customers_copy[r].demand &&
            vehicles_copy[j].time + t <= customers_copy[r].time_window_end &&
            vehicles_copy[j].time + t + customers_copy[r].service_time < ret_to_dep)
        {
            vehicles_copy[j].capacity -= customers_copy[r].demand;
            vehicles_copy[j].time += t + customers_copy[r].service_time;
            vehicles_copy[j].x_cord = customers_copy[r].x_cord;
            vehicles_copy[j].y_cord = customers_copy[r].y_cord;
            route.push_back(customers_copy[r].id);
            customers_copy.erase(customers_copy.begin() + i);
            coefficients.clear();
        }
        else
        {
            if (i == customers_copy.size())
            {
                i = 0;
                result.push_back(route);
                if (++j == vehicles_copy.size() && customers_copy.size() != 0)
                {
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

int grasp(Transport transport, Customers customers)
{
    int best = -1;
    bool stop = false;
    while (!stop)
    {
        vector<vector<int>> instance = greedy_randomized(transport, customers);
        while (instance.empty())
        {
            instance = greedy_randomized(transport, customers);
        }
        int local_best = local_search(instance, transport, customers);
        if (local_best < best || best == -1)
        {
            best = local_best;
        }
    }
    return best;
}

int local_search(vector<vector<int>> solution, Transport transport, Customers customers)
{
    vector<vector<int>> new_solution = solution;
    int best_cost = cost(solution, transport, customers);
    while (1)
    {
        int i = rand() % new_solution.size();
        int j = rand() % new_solution.size();
        while (i == j)
        {
            j = rand() % new_solution.size();
        }

        int k = rand() % new_solution[i].size();
        int l = rand() % new_solution[j].size();

        swap(new_solution[i][k], new_solution[j][l]);
        int new_cost = cost(new_solution, transport, customers);
        if (new_cost != -1 && new_cost < best_cost)
        {
            best_cost = new_cost;
        }
        else
        {
            swap(new_solution[i][k], new_solution[j][l]);
        }
    }
    return best_cost;
}

int cost(vector<vector<int>> solution, Transport transport, Customers customers)
{
    int total_cost = 0;
    for (int i = 0; i < solution.size(); i++)
    {
        vehicle dummy_vehicle(transport.vehicle_cap, transport.x_cord_of_dispatcher, transport.y_cord_of_dispatcher);
        for (int j = 0; j < solution[i].size(); j++)
        {
            int t = time_to_arrive(dummy_vehicle, customers.customers[i]);
            if (t > customers.customers[i].time_window_end)
            {
                return -1;
            }
            dummy_vehicle.capacity -= customers.customers[i].demand;
            dummy_vehicle.time += t + customers.customers[i].service_time;
            dummy_vehicle.x_cord = customers.customers[i].x_cord;
            dummy_vehicle.y_cord = customers.customers[i].y_cord;
        }
        int t = time_to_arrive(dummy_vehicle, customers.customers[0]);
        if (t > customers.customers[0].time_window_end)
        {
            return -1;
        }
        total_cost += dummy_vehicle.time + t;
    }
    return total_cost;
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
