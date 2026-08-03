#pragma once

#include <material.hpp>

class Layer {
    Material _material;
    double _thickness;

    public:
        Layer(Material& material, double thickness);

        const Material& getMaterial() const;
        const double getThickness() const;
};