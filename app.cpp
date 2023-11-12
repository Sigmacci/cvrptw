#include <iostream>
#include <vector>
#include <cmath>

#define VEHICLE_PARAMS 4
#define RECEIVER_PARAMS 7

#define VEHICLE_CAPACITY_PARAM 0
#define VEHICLE_TIME_PARAM 1
#define VEHICLE_X_CORD_PARAM 2
#define VEHICLE_Y_CORD_PARAM 3

#define RECEIVER_ID 0
#define RECEIVER_X_CORD_PARAM 1
#define RECEIVER_Y_CORD_PARAM 2
#define RECEIVER_DEMAND_PARAM 3
#define RECEIVER_TIME_WINDOW_START_PARAM 4
#define RECEIVER_TIME_WINDOW_END_PARAM 5
#define RECEIVER_SERVICE_TIME_PARAM 6

int time_to_arrive(std::vector<int> vehicle, std::vector<int> receiver);
void heapify(std::vector<std::vector<int>> M, std::vector<double> coefficient, int n, int i);
void heapsort(std::vector<std::vector<int>> M, std::vector<double> coefficient, int n);
int greedy(std::vector<std::vector<int>> vehicles, std::vector<std::vector<int>> receivers);

int main(int argc, char *argv[])
{
    std::vector<std::vector<int>> vehicles;
    std::vector<std::vector<int>> receivers;
}

int time_to_arrive(std::vector<int> vehicle, std::vector<int> receiver)
{
    return static_cast<int>(sqrt(pow(vehicle[VEHICLE_X_CORD_PARAM] - receiver[RECEIVER_X_CORD_PARAM], 2) +
            pow(vehicle[VEHICLE_Y_CORD_PARAM] - receiver[RECEIVER_Y_CORD_PARAM], 2)));
}

void heapify(std::vector<std::vector<int>> M, std::vector<double> coefficient, int n, int i)
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
		std::swap(coefficient[i], coefficient[min]);
		std::swap(M[i], M[min]);
		heapify(M, coefficient, n, min);
	}
}

void heapsort(std::vector<std::vector<int>> M, std::vector<double> coefficient, int n)
{
	for (int i = n / 2 - 1; i >= 0; i--)
	{
		heapify(M, coefficient, n, i);
	}
	for (int i = 1; i < n; i++)
	{
		std::swap(coefficient[0], coefficient[n - i]);
		std::swap(M[0], M[n - i]);
		heapify(M, coefficient, n - i, 0);
	}
}

int greedy(std::vector<std::vector<int>> vehicles, std::vector<std::vector<int>> receivers)
{
    std::vector<double> coefficients;
    std::vector<std::vector<int>> receivers_copy;
    std::copy(receivers, receivers_copy, receivers.size());
    int i = 0, j = 0;
    while (i < receivers_copy.size())
    {
        for (int k = 0; k < receivers_copy.size(); k++)
        {
            coefficients.push_back(time_to_arrive(vehicles[j], receivers_copy[k]) / receivers_copy[k][RECEIVER_DEMAND_PARAM]);
        }
        heapsort(receivers_copy, coefficients, receivers_copy.size());
        if (vehicles[j][VEHICLE_TIME_PARAM] < receivers_copy[i][RECEIVER_TIME_WINDOW_START_PARAM])
        {
            vehicles[j][VEHICLE_TIME_PARAM] = receivers_copy[i][RECEIVER_TIME_WINDOW_START_PARAM];
        }
        else if (vehicles[j][VEHICLE_TIME_PARAM] > receivers_copy[i][RECEIVER_TIME_WINDOW_END_PARAM])
        {
            j++;
            if (j == vehicles.size())
            {
                return 0;
            }
            continue;
        }
        if (vehicles[j][VEHICLE_CAPACITY_PARAM] >= receivers_copy[i][RECEIVER_DEMAND_PARAM] &&
            vehicles[j][VEHICLE_TIME_PARAM] + time_to_arrive(vehicles[j], receivers_copy[i]) <= receivers_copy[i][RECEIVER_TIME_WINDOW_END_PARAM])
        {
            vehicles[j][VEHICLE_CAPACITY_PARAM] -= receivers_copy[i][RECEIVER_DEMAND_PARAM];
            vehicles[j][VEHICLE_TIME_PARAM] += time_to_arrive(vehicles[j], receivers_copy[i]) + receivers_copy[i][RECEIVER_SERVICE_TIME_PARAM];
            vehicles[j][VEHICLE_X_CORD_PARAM] = receivers_copy[i][RECEIVER_X_CORD_PARAM];
            vehicles[j][VEHICLE_Y_CORD_PARAM] = receivers_copy[i][RECEIVER_Y_CORD_PARAM];
            receivers_copy.erase(receivers_copy.begin() + i);
            coefficients.clear();
            i++;
        }
    }
    int result = 0;
    for (auto vehicle : vehicles)
    {
        result += vehicle[VEHICLE_TIME_PARAM];
    }
    return result;
}