#include "controlview.h"

using namespace univiewer;

#if defined(WIN321) && !defined(_DEBUG) 
int WinMain(HINSTANCE, HINSTANCE, LPSTR, int) {

  std::vector<std::string> argv_str;
  int argc;
  LPWSTR *argv = CommandLineToArgvW(GetCommandLineW(), &argc);
  if (argc) {
    for (int i = 1; i < argc; ++i) {
      char str[1024];
      size_t retval;
      wcstombs_s(&retval, str, sizeof(str), argv[i], sizeof(str));
      argv_str.push_back(str);
    }
    LocalFree(argv);
  }
#else
int main(int argc, char *argv[]) {
  std::vector<std::string> argv_str;

  for (int i = 1; i < argc; ++i) {
    argv_str.push_back(argv[i]);
  }

  // argv_str.push_back("/5");
  // argv_str.push_back("D:/repos/omega-py/data/sln_rigidbody3.hdf5");

#endif

  sptr<ControlView> pControlView = CreateOneOf<ControlView>();
  pControlView->SetRender();

  if (!argv_str.empty())
    pControlView->InputModelfiles(argv_str);

  pControlView->SetAnimationMethod(DEFAULT_TIMERCALLBACK);
  pControlView->SetKeyboardMethod(DEFAULT_KEYPRESSCALLBACK);
  pControlView->SetWindowMethod(DEFAULT_WINDOWCALLBACK);
  pControlView->Display();

  return EXIT_SUCCESS;
}

