#include "es/config.hpp"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "es/constants.hpp"

namespace es {
namespace {

std::string trim(const std::string& s) {
  std::size_t a = 0, b = s.size();
  while (a < b && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
  while (b > a && std::isspace(static_cast<unsigned char>(s[b - 1]))) --b;
  return s.substr(a, b - a);
}

std::string lower(std::string s) {
  for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  return s;
}

struct Unit {
  const char* suffix;
  Real factor;
  Real offset;
};

// Ordered longest-first so that "mm" is not eaten by "m".
const Unit kUnits[] = {
    {"nL/s", 1e-12, 0}, {"uL/s", 1e-9, 0},  {"pL/s", 1e-15, 0}, {"mL/s", 1e-6, 0},
    {"m3/s", 1.0, 0},
    {"mbar", 100.0, 0}, {"kPa", 1e3, 0},    {"bar", 1e5, 0},    {"Pa", 1.0, 0},
    {"deg", constants::pi / 180.0, 0},      {"rad", 1.0, 0},
    {"kV", 1e3, 0},
    {"pA", 1e-12, 0},   {"nA", 1e-9, 0},    {"uA", 1e-6, 0},    {"mA", 1e-3, 0},
    {"nm", 1e-9, 0},    {"um", 1e-6, 0},    {"mm", 1e-3, 0},    {"cm", 1e-2, 0},
    {"S/m", 1.0, 0},    {"N/m", 1.0, 0},
    {"A", 1.0, 0},      {"V", 1.0, 0},      {"K", 1.0, 0},      {"C", 1.0, 273.15},
    {"m", 1.0, 0},
};

}  // namespace

Real Config::parse_value(const std::string& text) {
  const std::string t = trim(text);
  if (t.empty()) throw std::runtime_error("empty numeric value");

  // Split into the numeric head and the unit tail.
  std::size_t i = 0;
  if (t[i] == '+' || t[i] == '-') ++i;
  while (i < t.size() && (std::isdigit(static_cast<unsigned char>(t[i])) || t[i] == '.')) ++i;
  if (i < t.size() && (t[i] == 'e' || t[i] == 'E')) {
    // Exponent, but only if it is really one: e12, e-3.  "e" alone is a unit-ish
    // typo and must not swallow the rest.
    std::size_t j = i + 1;
    if (j < t.size() && (t[j] == '+' || t[j] == '-')) ++j;
    if (j < t.size() && std::isdigit(static_cast<unsigned char>(t[j]))) {
      i = j;
      while (i < t.size() && std::isdigit(static_cast<unsigned char>(t[i]))) ++i;
    }
  }
  const std::string head = t.substr(0, i);
  const std::string tail = trim(t.substr(i));
  if (head.empty()) throw std::runtime_error("cannot parse number from '" + t + "'");

  const Real x = std::strtod(head.c_str(), nullptr);
  if (tail.empty()) return x;

  for (const Unit& u : kUnits)
    if (tail == u.suffix) return x * u.factor + u.offset;

  throw std::runtime_error("unknown unit '" + tail + "' in value '" + t + "'");
}

void Config::set(const std::string& key, const std::string& value) {
  kv_[lower(trim(key))] = trim(value);
}

void Config::load(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open config file " + path);
  std::string line;
  std::string section;
  int lineno = 0;
  while (std::getline(in, line)) {
    ++lineno;
    const std::size_t hash = line.find_first_of("#;");
    if (hash != std::string::npos) line = line.substr(0, hash);
    line = trim(line);
    if (line.empty()) continue;
    if (line.front() == '[' && line.back() == ']') {
      section = lower(trim(line.substr(1, line.size() - 2)));
      continue;
    }
    const std::size_t eq = line.find('=');
    if (eq == std::string::npos)
      throw std::runtime_error(path + ":" + std::to_string(lineno) + ": expected key = value");
    std::string key = lower(trim(line.substr(0, eq)));
    if (!section.empty() && key.find('.') == std::string::npos) key = section + "." + key;
    kv_[key] = trim(line.substr(eq + 1));
  }
}

std::vector<std::string> Config::positional_args(int argc, char** argv) {
  std::vector<std::string> rest;
  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    const std::size_t eq = a.find('=');
    if (eq == std::string::npos || eq == 0) rest.push_back(a);
  }
  return rest;
}

void Config::apply_cli(int argc, char** argv) {
  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    const std::size_t eq = a.find('=');
    if (eq == std::string::npos || eq == 0) continue;
    set(a.substr(0, eq), a.substr(eq + 1));
  }
}

bool Config::has(const std::string& key) const { return kv_.count(lower(key)) > 0; }

Real Config::num(const std::string& key, Real fallback) const {
  const std::string k = lower(key);
  read_.insert(k);
  const auto it = kv_.find(k);
  if (it == kv_.end()) return fallback;
  try {
    return parse_value(it->second);
  } catch (const std::exception& e) {
    throw std::runtime_error("config key '" + k + "': " + e.what());
  }
}

int Config::integer(const std::string& key, int fallback) const {
  return static_cast<int>(std::lround(num(key, static_cast<Real>(fallback))));
}

bool Config::flag(const std::string& key, bool fallback) const {
  const std::string k = lower(key);
  read_.insert(k);
  const auto it = kv_.find(k);
  if (it == kv_.end()) return fallback;
  const std::string v = lower(it->second);
  if (v == "true" || v == "yes" || v == "on" || v == "1") return true;
  if (v == "false" || v == "no" || v == "off" || v == "0") return false;
  throw std::runtime_error("config key '" + k + "': expected a boolean, got '" + it->second + "'");
}

std::string Config::str(const std::string& key, const std::string& fallback) const {
  const std::string k = lower(key);
  read_.insert(k);
  const auto it = kv_.find(k);
  return (it == kv_.end()) ? fallback : it->second;
}

std::vector<std::string> Config::unused_keys(const std::vector<std::string>& ignore) const {
  std::vector<std::string> v;
  for (const auto& kv : kv_) {
    if (read_.count(kv.first)) continue;
    bool skip = false;
    for (const std::string& pre : ignore)
      if (kv.first.rfind(lower(pre), 0) == 0) { skip = true; break; }
    if (!skip) v.push_back(kv.first);
  }
  return v;
}

void Config::warn_about_unused(std::FILE* out, const std::vector<std::string>& ignore) const {
  const std::vector<std::string> u = unused_keys(ignore);
  if (u.empty()) return;
  std::fprintf(out, "\nWARNING: %zu configuration key(s) were never used.  A misspelled key is\n"
                    "         silently ignored otherwise, so check these:\n", u.size());
  for (const std::string& k : u) std::fprintf(out, "         %s = %s\n", k.c_str(),
                                              kv_.at(k).c_str());
}

void Config::dump(std::FILE* out) const {
  for (const auto& kv : kv_) std::fprintf(out, "  %-34s = %s\n", kv.first.c_str(), kv.second.c_str());
}

}  // namespace es
