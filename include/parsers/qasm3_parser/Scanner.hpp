#pragma once

#include "Token.hpp"

#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stack>

namespace qasm3 {
class Scanner {
private:
  std::istream& is;
  std::map<std::string, Token::Kind, std::less<>> keywords{};
  char ch = 0;
  size_t line = 1;
  size_t col = 0;

  [[nodiscard]] static bool isSpace(const char c) {
    return c == ' ' || c == '\t' || c == '\r' || c == '\n';
  }

  [[nodiscard]] static bool isFirstIdChar(const char c) {
    return isalpha(c) != 0 || c == '_';
  }

  [[nodiscard]] static bool isNum(const char c) { return c >= '0' && c <= '9'; }

  [[nodiscard]] static bool isHex(const char c) {
    return isNum(c) || (c >= 'a' && c <= 'f') || (c >= 'A' && c <= 'F');
  }

  static char readUtf8Codepoint(std::istream& ss);

  void nextCh();

  [[nodiscard]] char peek() const;

  std::optional<Token> consumeWhitespaceAndComments();

  uint64_t consumeIntegerLiteral(uint8_t base);

  Token consumeNumberLiteral();

  Token consumeHardwareQubit();

  Token consumeString();

  Token consumeName();

  void error(const std::string& msg) const {
    std::cerr << "Error at line " << line << ", column " << col << ": " << msg
              << '\n';
  }

  void expect(const char expected) {
    if (ch != expected) {
      error("Expected '" + std::to_string(expected) + "', got '" + ch + "'");
    } else {
      nextCh();
    }
  }

public:
  explicit Scanner(std::istream& in);

  ~Scanner() = default;

  Token next();
};
} // namespace qasm3