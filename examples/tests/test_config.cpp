#include <configuration.hpp>

#include <iostream>
#include <matplot/matplot.h>

int main() {
    SimulationManager manager("/src/examples/tests/simple_emitter.json");
    auto solver = manager.create();

    solver->calculate();


    // Polar figure
    Eigen::ArrayXd thetaGlass, powerPerpAngleGlass, powerParasPolAngleGlass, powerParapPolAngleGlass;
    solver->calculateEmissionSubstrate(thetaGlass, powerPerpAngleGlass, powerParapPolAngleGlass, powerParasPolAngleGlass);

    std::cout << solver->mPowerPerpUpPol.leftCols(5) << '\n';

    matplot::figure();
    matplot::plot(thetaGlass, powerPerpAngleGlass, "-o");
    matplot::hold(matplot::on);
    matplot::plot(thetaGlass, powerParapPolAngleGlass, "-o");
    matplot::plot(thetaGlass, powerParasPolAngleGlass, "-o");
    
    matplot::figure();
    matplot::polarplot(thetaGlass, powerPerpAngleGlass, "-")->line_width(2);
    matplot::hold(matplot::on);
    matplot::polarplot(thetaGlass, powerParapPolAngleGlass, "-")->line_width(2);
    matplot::polarplot(thetaGlass, powerParasPolAngleGlass, "-")->line_width(2);

    matplot::show();
}