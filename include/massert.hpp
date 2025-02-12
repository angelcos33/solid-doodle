// Simple macro for custom assert with message
#pragma once

#include <cassert>

#define m_assert(expr, msg) assert(((void)(msg), (expr)))