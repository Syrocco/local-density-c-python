#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <omp.h>
#include <random>

namespace py = pybind11;

class CellList {
private:
    std::vector<double> x, y, r;
    double Lx, Ly, cell_size;
    std::vector<std::vector<int>> cells;
    int nx, ny;

public:
    int toCellX(double x_val) {
        return static_cast<int>(x_val / cell_size)%nx;
    }

    int toCellY(double y_val) {
        return static_cast<int>(y_val / cell_size)%ny ;
    }

    double pbcX(double dx) {
        if (dx >= Lx / 2)
            return dx - Lx;
        else if (dx <= -Lx / 2)
            return dx + Lx;
        return dx;
    }

    double pbcY(double dy) {
        if (dy >= Ly / 2)
            return dy - Ly;
        else if (dy <= -Ly / 2)
            return dy + Ly;
        return dy;
    }

    int pbcCellX(int a) {
        if (a < 0)
            return a + nx;
        else if (a >= nx)
            return a - nx;
        return a;
    }

    int pbcCellY(int a) {
        if (a < 0)
            return a + ny;
        else if (a >= ny)
            return a - ny;
        return a;
    }

    CellList(std::vector<double> x_vals, std::vector<double> y_vals, std::vector<double> r_vals, double Lx_val, double Ly_val, double cell_size_val) {
        x = x_vals;
        y = y_vals;
        Lx = Lx_val;
        Ly = Ly_val;
        r = r_vals;
        cell_size = cell_size_val;
        nx = static_cast<int>(Lx / cell_size);
        ny = static_cast<int>(Ly / cell_size);
        cells.resize(nx * ny);

        for (size_t idx = 0; idx < x.size(); ++idx) {
            int cell_x = toCellX(x[idx]);
            int cell_y = toCellY(y[idx]);
            cells[cell_y*nx + cell_x].push_back(idx);
        }
    }

    std::vector<int> get_particles_in_cell(int cell_x, int cell_y) {
        return cells[cell_y*nx + cell_x];
    }

    double get_density(double x_val, double y_val, double rad, double radius = 0.56123102415) {
        double mini = pow(rad - radius, 2);
        double maxi = pow(rad + radius, 2);
        auto f = [mini, maxi, rad, radius](double X) {
            if (X < mini)
                return 1.0;
            else if (X > maxi)
                return 0.0;
            else{
                double xx = sqrt(X);
                double angle1 = acos((radius*radius + X - rad*rad)/(2*radius*xx));
                double angle2 = acos((rad*rad + X - radius*radius)/(2*rad*xx));
                double area1 = radius*radius*angle1 - radius*radius*sin(2*angle1)/2;
                double area2 = rad*rad*angle2 - rad*rad*sin(2*angle2)/2;
                return (area1 + area2)/(M_PI*radius*radius);
                
            }
        };

        int numCell = static_cast<int>(rad / cell_size) + 2;
        int cell_x_i = toCellX(x_val);
        int cell_y_i = toCellY(y_val);
        double c = 0.0;
        for (int i = -numCell; i < numCell; ++i) {
            for (int j = -numCell; j < numCell; ++j) {
                auto particles_in_cell = get_particles_in_cell(pbcCellX(cell_x_i + i), pbcCellY(cell_y_i + j));
                for (auto k : particles_in_cell) {
                    double distance_squared = pow(pbcX(x_val - x[k]), 2) + pow(pbcY(y_val - y[k]), 2);
                    c += f(distance_squared);
                }
            }
        }
        return c / (M_PI * pow(rad, 2))*M_PI*pow(radius, 2) ;
    }


    double get_rad_density(double x_val, double y_val, double rad) {
        
        auto f = [rad](double X, double radius) {
            double mini = pow(rad - radius, 2);
            double maxi = pow(rad + radius, 2);
            if (X < mini)
                return M_PI*radius*radius;
            else if (X > maxi)
                return 0.0;
            else{
                double xx = sqrt(X);
                double angle1 = acos((radius*radius + X - rad*rad)/(2*radius*xx));
                double angle2 = acos((rad*rad + X - radius*radius)/(2*rad*xx));
                double area1 = radius*radius*angle1 - radius*radius*sin(2*angle1)/2;
                double area2 = rad*rad*angle2 - rad*rad*sin(2*angle2)/2;
                return area1 + area2;
                
            }
        };

        int numCell = static_cast<int>(rad / cell_size) + 2;
        int cell_x_i = toCellX(x_val);
        int cell_y_i = toCellY(y_val);
        double c = 0.0;
        for (int i = -numCell; i < numCell; ++i) {
            for (int j = -numCell; j < numCell; ++j) {
                auto particles_in_cell = get_particles_in_cell(pbcCellX(cell_x_i + i), pbcCellY(cell_y_i + j));
                for (auto k : particles_in_cell) {
                    double distance_squared = pow(pbcX(x_val - x[k]), 2) + pow(pbcY(y_val - y[k]), 2);
               
                    c += f(distance_squared, r[k]);
                }
            }
        }
        return c / (M_PI * pow(rad, 2)) ;
    }

    double get_linear_density(double x_val, double y_val, double rad, double radius = 0.56123102415) {
        double mini = pow(rad - radius, 2);
        double maxi = pow(rad + radius, 2);
        auto f = [mini, maxi, rad, radius](double X) {
            if (X < mini)
                return 1.0;
            else if (X > maxi)
                return 0.0;
            else{
                return (1.0 - 0.0) / (mini - maxi) * (X - maxi);
            }
        };

        int numCell = static_cast<int>(rad / cell_size) + 2;
        int cell_x_i = toCellX(x_val);
        int cell_y_i = toCellY(y_val);
        double c = 0.0;
        for (int i = -numCell; i < numCell; ++i) {
            for (int j = -numCell; j < numCell; ++j) {
                auto particles_in_cell = get_particles_in_cell(pbcCellX(cell_x_i + i), pbcCellY(cell_y_i + j));
                for (auto k : particles_in_cell) {
                    double distance_squared = pow(pbcX(x_val - x[k]), 2) + pow(pbcY(y_val - y[k]), 2);
                    c += f(distance_squared);
                }
            }
        }
        return c / (M_PI * pow(rad, 2))*M_PI*pow(radius, 2) ;
    }

    double get_average_density(int N, double rad, double radius = 0.56123102415){
        std::random_device r;
        std::vector<std::default_random_engine> generators;
        for (int i = 0, N = omp_get_max_threads(); i < N; ++i) {
            generators.emplace_back(std::default_random_engine(r()));
        }
        std::uniform_real_distribution<float> distributionX(0.0, Lx);
        std::uniform_real_distribution<float> distributionY(0.0, Ly);


        double average = 0;
        #pragma omp parallel for
        for (int i = 0; i < N; i++){
            std::default_random_engine& engine = generators[omp_get_thread_num()];
            
            double x_rand = distributionX(engine);
            double y_rand = distributionY(engine);
            average += get_density(x_rand, y_rand, rad, radius);
        }
        return average/N;
    }

    std::vector<double> get_density_list(int N, double rad,  double radius = 0.56123102415){
        std::random_device r;
        std::vector<std::default_random_engine> generators;
        for (int i = 0, N = omp_get_max_threads(); i < N; ++i) {
            generators.emplace_back(std::default_random_engine(r()));
        }
        std::uniform_real_distribution<float> distributionX(0.0, Lx);
        std::uniform_real_distribution<float> distributionY(0.0, Ly);

        std::vector<double> list(N);
        #pragma omp parallel for
        for (int i = 0; i < N; i++){
            std::default_random_engine& engine = generators[omp_get_thread_num()];
            
            double x_rand = distributionX(engine);
            double y_rand = distributionY(engine);
            list[i] = get_density(x_rand, y_rand, rad, radius);
        }
        return list;
    }

    std::vector<double> get_rad_density_list(int N, double rad){
        std::random_device r;
        std::vector<std::default_random_engine> generators;
        for (int i = 0, N = omp_get_max_threads(); i < N; ++i) {
            generators.emplace_back(std::default_random_engine(r()));
        }
        std::uniform_real_distribution<float> distributionX(0.0, Lx);
        std::uniform_real_distribution<float> distributionY(0.0, Ly);

        std::vector<double> list(N);
        #pragma omp parallel for
        for (int i = 0; i < N; i++){
            std::default_random_engine& engine = generators[omp_get_thread_num()];
            
            double x_rand = distributionX(engine);
            double y_rand = distributionY(engine);
            list[i] = get_rad_density(x_rand, y_rand, rad);
        }
        return list;
    }

    std::vector<double> get_linear_density_list(int N, double rad,  double radius = 0.56123102415){
        std::random_device r;
        std::vector<std::default_random_engine> generators;
        for (int i = 0, N = omp_get_max_threads(); i < N; ++i) {
            generators.emplace_back(std::default_random_engine(r()));
        }
        std::uniform_real_distribution<float> distributionX(0.0, Lx);
        std::uniform_real_distribution<float> distributionY(0.0, Ly);

        std::vector<double> list(N);
        #pragma omp parallel for
        for (int i = 0; i < N; i++){
            std::default_random_engine& engine = generators[omp_get_thread_num()];
            
            double x_rand = distributionX(engine);
            double y_rand = distributionY(engine);
            list[i] = get_linear_density(x_rand, y_rand, rad, radius);
        }
        return list;
    }
    py::array_t<double> get_linear_density_array(int N, double rad, double radius = 0.56123102415) {
        auto result = get_linear_density_list(N, rad, radius);
        py::array_t<double> arr(result.size());
        auto buffer = arr.request();
        double *ptr = static_cast<double *>(buffer.ptr);
        std::copy(result.begin(), result.end(), ptr);
        return arr;
    }
    py::array_t<double> get_density_array(int N, double rad, double radius = 0.56123102415) {
        auto result = get_density_list(N, rad, radius);
        py::array_t<double> arr(result.size());
        auto buffer = arr.request();
        double *ptr = static_cast<double *>(buffer.ptr);
        std::copy(result.begin(), result.end(), ptr);
        return arr;
    }

    py::array_t<double> get_rad_density_array(int N, double rad) {
        auto result = get_rad_density_list(N, rad);
        py::array_t<double> arr(result.size());
        auto buffer = arr.request();
        double *ptr = static_cast<double *>(buffer.ptr);
        std::copy(result.begin(), result.end(), ptr);
        return arr;
    }
};




PYBIND11_MODULE(cell_list, m) {
    py::class_<CellList>(m, "CellList")
        .def(py::init<std::vector<double>, std::vector<double>, std::vector<double>, double, double, double>())
        .def("get_density", &CellList::get_density)
        .def("get_average_density", &CellList::get_average_density)
        .def("get_linear_density_array", &CellList::get_linear_density_array)
        .def("get_density_array", &CellList::get_density_array)
        .def("get_rad_density_array", &CellList::get_rad_density_array);
}