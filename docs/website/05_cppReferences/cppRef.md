<!--
SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
SPDX-License-Identifier: CC-BY-SA-4.0
-->

# C++ recommendations
Since Ikarus is written in C++, we summarize our recommendations to dig deeper into C++ coding on this page.

## Should I write a member function or a free function?
1. Scott Meyers recommends the algorithm in the [link](http://cpptips.com/nmemfunc_encap).

## When to use const?
1. Refer the blog of [Arthur O'Dwyer](https://quuxplusone.github.io/blog/2022/01/23/dont-const-all-the-things/)

## How should I pass my parameters to a function and return from a function?
1. Refer [Herb Sutter's Cppcon Talk 2014](https://youtu.be/xnqTKD8uD64?t=3318)

## Best practices

1. [C++ Core Guidelines](https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines)
2. [Jason Turner's collection of best practices](https://lefticus.gitbooks.io/cpp-best-practices)
3. [More C++ idioms](https://en.wikibooks.org/wiki/More_C%2B%2B_Idioms)

## Further references

1. [Cppcon Videos](https://www.youtube.com/user/CppCon) - These videos are released after every C++ conference. For beginners, the "Back to basics" track is recommended.
2. [Godblot](https://godbolt.org/) - Online compiler with assembler output. It's useful to quickly determine whether something will be fast or slow.
   Furthermore, libraries like Eigen can be added. Also, any other header files found on the internet can be included with the link.

## Videos
Here we collect some useful videos on general coding or coding with C++:

1. [Clean Code - Uncle Bob / Lesson 1](https://www.youtube.com/watch?v=7EmboKQH8lM) How to write code cleanly, see also [@martinclean]
2. [CppCon 2014: Herb Sutter "Back to the Basics! Essentials of Modern C++ Style"](https://youtu.be/xnqTKD8uD64)
3. [CppCon 2018: Jonathan Boccara “105 STL Algorithms in Less Than an Hour”](https://www.youtube.com/watch?v=2olsGf6JIkU) "Almost" all algorithms in the STL
4. [Back to Basics: Object-Oriented Programming - Jon Kalb - CppCon 2019](https://www.youtube.com/watch?v=32tDTD9UJCE&list=PLHTh1InhhwT4CTnVjJqnAKeMfGzOWjsRa) How to do modern "Object-Oriented Programming" (If you really have to)
5. [CppCon 2021 - Back To Basics](https://www.youtube.com/watch?v=Bt3zcJZIalk&list=PLHTh1InhhwT4TJaHBVWzvBOYhp27UO7mI)
6. [CppCon 2019 Back to Basics](https://www.youtube.com/watch?v=32tDTD9UJCE&list=PLHTh1InhhwT4CTnVjJqnAKeMfGzOWjsRa)

## Books
1. [Meyers S. 1995][@meyers1995more]
2. [Gamma E. et al. 1995][@gamma1995design]
3. [Meyers S. 2005][@meyers2005effective]
4. [Reddy M. 2011][@reddy2011api]
5. [Iglberger K. 2022][@iglbergerSoftwareDesignDesign2022]

\bibliography 