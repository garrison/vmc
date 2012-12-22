solution "vmc-core"
  configurations { "Local", "Production", "Debug", "Careful" }

  project "vmc-core"
    language "C++"
    kind "SharedLib"
    files { "*.hpp", "*.cpp" }

    -- you can add to these by passing INCLUDES as an environment
    -- variable.  see the README for an example.
    includedirs { "/usr/include/eigen3", "/usr/include/boost" }

    buildoptions { "-Wall", "-Wextra" }

    -- see http://shorestreet.com/why-your-dso-is-slow
    configuration "linux"
      linkoptions { "-Wl,-Bsymbolic-functions" }

    configuration "Local"
      flags { "OptimizeSpeed" }

    configuration "Production"
      defines { "NDEBUG" }
      flags { "OptimizeSpeed" }

    configuration "Debug"
      flags { "Symbols" }
      flags { "Optimize" }

    configuration "Careful"
      flags { "Symbols" }
      defines { "VMC_CAREFUL" }
      flags { "Optimize" }
