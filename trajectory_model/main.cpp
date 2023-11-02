#include "trajectory.h"

int main()
{
    PredictPitchXY &trajectory = PredictPitchXY::getinstance();
    trajectory(16, 20, 2);
    return 0;
}