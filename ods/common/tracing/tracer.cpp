#include "common/tracing/tracer.hpp"

#ifndef ENCLAVE_MODE

static void print_stacktrace(int signum) {
  ::signal(signum, SIG_DFL);
  std::cerr << "Stack trace:\n" << boost::stacktrace::stacktrace() << '\n';
  Profiler::g_Tracker.Reset();
  TimeTracker::g_Tracker.Reset();
  g_MaxTracker.Reset();
  ::raise(SIGABRT);
}

struct OnExitHandlerInstaller {
  int theStart;
  OnExitHandlerInstaller() {
    theStart = 10;
    std::cerr << "Installing signal handler" << std::endl;
    ::signal(SIGSEGV, &print_stacktrace);
    ::signal(SIGABRT, &print_stacktrace);
  }
};

TimeTracker TimeTracker::g_Tracker;
Profiler Profiler::g_Tracker;

MaxTracker g_MaxTracker;
bool g_disableTracing = false;
bool g_disableProfiling = true;
OnExitHandlerInstaller g_OnExit;

#endif