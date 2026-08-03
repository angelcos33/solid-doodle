#include <utility>

#include <layer.hpp>

Layer::Layer(Material& material, double thickness):
             _material(std::move(material)),
             _thickness(thickness)
{}

const Material& Layer::getMaterial() const {
    return _material;
}

const double Layer::getThickness() const {
    return _thickness;
}