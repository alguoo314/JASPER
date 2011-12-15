#include <iostream>
#include <test.hpp>
#include <src/test_cmdline.hpp>

int main(int argc, char *argv[]) {
  test_cmdline args(argc, argv);
  std::cout << "Hello the world!\n";
  for(test_cmdline::name_arg_it it = args.name_arg.begin(); it != args.name_arg.end(); ++it)
    std::cout << "Hello " << *it << "\n";

  return 0;
}
