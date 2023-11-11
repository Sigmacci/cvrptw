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

int eval_time(std::vector<int> vehicle, std::vector<int> receiver)
{
    return vehicle[VEHICLE_TIME_PARAM] + 
        sqrt(pow(vehicle[VEHICLE_X_CORD_PARAM] - receiver[RECEIVER_X_CORD_PARAM], 2) +
            pow(vehicle[VEHICLE_Y_CORD_PARAM] - receiver[RECEIVER_Y_CORD_PARAM], 2));
}

void heapify(std::vector<std::vector<int>> M, double coefficient[], int n, int i)
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

void heapsort(std::vector<std::vector<int>> M, double coefficient[], int n)
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

int *greedy(std::vector<std::vector<int>> vehicles, std::vector<std::vector<int>> receivers)
{
    double coefficient[receivers.size()];
    for (int i = 0; i < receivers.size(); i++)
    {
        coefficient[i] = receivers[i][RECEIVER_DEMAND_PARAM] / eval_time(vehicles[0], receivers[i]);
    }
    heapsort(receivers, coefficient, receivers.size());
    int i = 0;
    int j = 0;
    int *k = new int[vehicles.size()];
    while (i < receivers.size())
    {
        if (vehicles[j][VEHICLE_CAPACITY_PARAM] >= receivers[i][RECEIVER_DEMAND_PARAM] &&
            vehicles[j][VEHICLE_TIME_PARAM] + eval_time(vehicles[j], receivers[i]) <= receivers[i][RECEIVER_TIME_WINDOW_END_PARAM])
        {
            vehicles[j][VEHICLE_CAPACITY_PARAM] -= receivers[i][RECEIVER_DEMAND_PARAM];
            vehicles[j][VEHICLE_TIME_PARAM] += eval_time(vehicles[j], receivers[i]) + receivers[i][RECEIVER_SERVICE_TIME_PARAM];
            i++;
        }
        else
        {
            j++;
            if (j == vehicles.size())
            {
                return k;
            }
        }
    }
}

int main(int argc, char *argv[])
{
    std::vector<std::vector<int>> vehicles;
    std::vector<std::vector<int>> receivers;
}