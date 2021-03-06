newoption {
  trigger = "gtest",
  description = "googletest source root directory (default /usr/src/gtest)",
  value = "path",
}

local gtest_root = _OPTIONS["gtest"] or "/usr/src/gtest"

solution "vmc-core-tests"
  configurations { "Local" }

  buildoptions { "-g", "-Wall", "-Wextra", "-std=c++11" }
  includedirs { gtest_root .. "/include" }
  defines { "GTEST_HAS_PTHREAD=0" }

  configuration "Local"
    flags { "Optimize" }

  project "gtest_main"
    language "C++"
    kind "StaticLib"

    includedirs { gtest_root }
    files { gtest_root .. "/src/gtest-all.cc", gtest_root .. "/src/gtest_main.cc" }

  project "vmc-core-tests"
    language "C++"
    kind "ConsoleApp"

    files { "*.cpp" }
    includedirs { "../vmc-core" }

    links { "vmc-core", "gtest_main" }
    libdirs { "../vmc-core" }
