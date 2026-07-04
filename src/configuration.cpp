#include <iostream>

#include <configuration.hpp>


ConfigurationManager::ConfigurationManager(const std::string& filename):
    jsonParser(filename)
    {}

void ConfigurationManager::configure() {
    // Parse Config file
    jsonParser.parse();

    // Materials
    auto jsonTree = jsonParser.getJsonTree();
    jsonTree->print();
    auto stack = jsonTree->find("stack");
    std::vector<Material> materials;
    std::vector<double> thicknesses;
    double emitterPosition;
    bool hasEmitter = false;
    Eigen::Index emitterIndex;
    if (stack.has_value()) {
        auto it = stack.value();
        if (std::holds_alternative<Json::JsonList*>(it->second->value)) {
            auto layerList = std::get<Json::JsonList*>(it->second->value);
            for (auto layerIt = layerList->begin(); layerIt != layerList->end(); ++layerIt) {
                if (std::holds_alternative<Json::JsonObject*>((*layerIt)->value)) {
                    auto layer = std::get<Json::JsonObject*>((*layerIt)->value);
                    auto layerObject = layer->begin();
                    auto layerData = std::get<Json::JsonObject*>(layerObject->second->value);
                    auto materialIt = layerData->find("material");
                    if (materialIt != layerData->end()) {
                        materials.emplace_back(std::get<double>(materialIt->second->value), 0.0);
                    }
                    else {throw std::runtime_error("Each layer must have a material!");}
                    if (layerObject->first != "Environment" && layerObject->first != "Substrate") {
                        auto thicknessIt = layerData->find("thickness");
                        if (thicknessIt != layerData->end()) {
                            if (std::holds_alternative<double>(thicknessIt->second->value)) {
                                thicknesses.push_back(std::get<double>(thicknessIt->second->value));
                            }
                            else {throw std::runtime_error("Thickness must be a numeric value!");}
                        }
                        else {throw std::runtime_error("Every internal layer must have a thickness!");}
                        if (layerObject->first == "Emitter") {
                            hasEmitter = true;
                            emitterIndex = layerIt - layerList->begin();
                            auto emitterIt = layerData->find("position"); 
                            if (emitterIt != layerData->end()) {
                               if (std::holds_alternative<double>(emitterIt->second->value)) {
                                emitterPosition = std::get<double>(emitterIt->second->value);
                            }
                                else {throw std::runtime_error("Emitter position must be a numeric value!");} 
                            }
                        }
                    }
                }
            }
        }
    }
    if (!hasEmitter) {throw std::runtime_error("You need at least one emitter!");}
    for (const auto& x: thicknesses) {
        std::cout << x << std::endl;
    }
    std::cout << "Emitter Position: " << emitterPosition << std::endl << "Emitter Index: " << emitterIndex << std::endl;
}
