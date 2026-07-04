#include <configuration.hpp>

int main() {
    ConfigurationManager manager("/src/examples/tests/simple_emitter.json");
    manager.configure();
}