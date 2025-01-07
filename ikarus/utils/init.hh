// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file init.hh
 * \brief Implementation of the init function
 */

#pragma once

#include <chrono>
#include <utility>

#include <dune/common/parallel/mpihelper.hh>

#include <spdlog/details/registry.h>
#include <spdlog/fmt/chrono.h>
#include <spdlog/fmt/ostr.h>
#include <spdlog/fmt/ranges.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

namespace Ikarus {

/**
 * \brief Singleton class representing an instance of the Ikarus framework.
 *
 * This class manages the initialization and configuration of the Ikarus framework.
 * It provides functionality to enable a file logger and prints a banner and version information on initialization.
 */
class IkarusInstance
{
public:
  /**
   * \brief Gets the singleton instance of the Ikarus framework.
   *
   * \return The singleton instance.
   */
  static IkarusInstance& getInstance() {
    static IkarusInstance instance;
    return instance;
  }

  /**
   * \brief Enables a file logger.
   *
   * \param filename The name of the log file. If empty, the executable name is used.
   */
  void enableFileLogger(std::string&& filename = "") {
    using namespace std::chrono;
    std::string currentTime = fmt::format("_{}", std::chrono::system_clock::now());

    std::ranges::transform(currentTime, currentTime.begin(), [](char ch) {
      return (ch == ' ' or ch == ':') ? '_' : ch;
    }); // replace space and colon with underscore
    auto logFilename = (filename.empty() ? executableName_ : filename) + currentTime + ".log";
    auto fileSink    = std::make_shared<spdlog::sinks::basic_file_sink_mt>(logFilename, true);
    fileSink->set_pattern("%v");
    fileSink->set_level(spdlog::level::trace);
    auto logger = spdlog::default_logger();
    logger->sinks().push_back(fileSink);
  }

private:
  friend void init(int argc, char** argv, bool enableFileLogger);
  IkarusInstance() = default;
  std::string executableName_;

public:
  IkarusInstance(const IkarusInstance&) = delete;
  void operator=(const IkarusInstance&) = delete;
};

/**
 * \brief Initializes the Ikarus framework.
 *
 * \param argc The number of command-line arguments.
 * \param argv The command-line arguments.
 * \param enableFileLogger Flag indicating whether to enable the file logger.
 */
void inline init(int argc, char** argv, bool enableFileLogger = true) {
  Dune::MPIHelper::instance(argc, argv);
  auto& instance           = IkarusInstance::getInstance();
  instance.executableName_ = argv[0];
  auto logger              = spdlog::default_logger();
  logger->set_pattern("%v");
  if (enableFileLogger)
    instance.enableFileLogger();

  auto currentTime = std::chrono::system_clock::now();
  const std::chrono::year_month_day currentYearMonthDay{floor<std::chrono::days>(currentTime)};

  spdlog::info("Start of execution: {}", currentTime);
  /// https://patorjk.com/software/taag/#p=testall&f=Univers&t=IKARUS (font: Lean)
  spdlog::info(R"xxx(    _/_/_/  _/    _/    _/_/    _/_/_/    _/    _/    _/_/_/)xxx");
  spdlog::info(R"xxx(     _/    _/  _/    _/    _/  _/    _/  _/    _/  _/       )xxx");
  spdlog::info(R"xxx(    _/    _/_/      _/_/_/_/  _/_/_/    _/    _/    _/_/    )xxx");
  spdlog::info(R"xxx(   _/    _/  _/    _/    _/  _/    _/  _/    _/        _/   )xxx");
  spdlog::info(R"xxx(_/_/_/  _/    _/  _/    _/  _/    _/    _/_/    _/_/_/      )xxx");

  spdlog::info("© 2021-{} The Ikarus Developers, see LICENSE.md ", static_cast<int>(currentYearMonthDay.year()));
  spdlog::info("You are using Ikarus v{}. Please don't forget to cite us:", IKARUS_VERSION);
  spdlog::info(
      "Müller, A., Vinod Kumar Mitruka, T. K. M., Jakob, H. (2024). Ikarus v0.4 (Version V1). "
      "doi:<https://doi.org/10.18419/darus-3889>");
}
} // namespace Ikarus
