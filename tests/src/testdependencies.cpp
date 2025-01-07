// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include <fstream>
#include <memory>

#include <dune/common/test/testsuite.hh>

#include <spdlog/spdlog.h>

#include <Eigen/Core>

#include <ikarus/utils/init.hh>

using Dune::TestSuite;

static void foo() {
  spdlog::info("Does this appear in the correct logger?");
  spdlog::debug("This is a debug statement");
}

static TestSuite spdlogTest() {
  TestSuite t("spdlogTest");

  auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
  console_sink->set_level(spdlog::level::trace);
  console_sink->set_pattern("[%^%l%$] %v");

  auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>("logs/testLog.txt", true);
  file_sink->set_level(spdlog::level::trace);
  file_sink->set_pattern("[%^%l%$] %v");

  // spdlog::logger logger("multi_sink", {console_sink, file_sink});
  spdlog::logger logger("multi_sink", {file_sink});
  spdlog::set_default_logger(std::make_shared<spdlog::logger>(logger));

  logger.info("Welcome to spdlog!");
  logger.error("Some error message with arg: {}", 1);

  logger.warn("Easy padding in numbers like {:08d}", 12);
  logger.critical("Support for int: {0:d};  hex: {0:x};  oct: {0:o}; bin: {0:b}", 42);
  logger.info("Support for floats {:03.2f}", 1.23456);
  logger.info("Positional args are {1} {0}..", "too", "supported");
  logger.info("{:<29}{}", "left aligned", "w");
  logger.info("{:>30}", "right aligned");

  logger.set_level(spdlog::level::debug); // Set global log level to debug
  logger.debug("This message should be displayed..");

  foo();
  Eigen::Vector3d a({1.0, 2.0, 3.0});
  logger.info("{}", a.transpose());
  logger.flush();
  std::ifstream ts(file_sink->filename());
  std::string file((std::istreambuf_iterator<char>(ts)), (std::istreambuf_iterator<char>()));
  std::string expectedOutput = R"xxx([info] Welcome to spdlog!
[error] Some error message with arg: 1
[warning] Easy padding in numbers like 00000012
[critical] Support for int: 42;  hex: 2a;  oct: 52; bin: 101010
[info] Support for floats 1.23
[info] Positional args are supported too..
[info] left aligned                 w
[info]                  right aligned
[debug] This message should be displayed..
[info] Does this appear in the correct logger?
[info] [1, 2, 3]
)xxx";
  t.check(file == expectedOutput);
  return t;
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t;

  t.subTest(spdlogTest());

  return t.exit();
}
