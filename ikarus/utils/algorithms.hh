// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file algorithms.hh
 * \brief Implementation of several stl-like algorithms

 */

#pragma once

// This file contains stl-like algorithms
#include <iosfwd>
#include <ranges>

#include <ikarus/utils/traits.hh>
namespace Ikarus::utils {

/**
 * \brief Sorts and removes duplicate elements from a random access range.
 * \ingroup algos
 *    *\details
 * This function sorts the elements of the given random access range and removes duplicate elements,
 * leaving only unique elements in the range.
 *
 * \param r The random access range to be modified.
 */
void makeUniqueAndSort(std::ranges::random_access_range auto& r) {
  sort(r.begin(), r.end());
  r.erase(std::unique(r.begin(), r.end()), r.end());
}

/**
 * \brief Appends a value to the range if it is not already present.
 *\details
 * This function appends a value to the given random access range only if the value is not already present in the
 *range. It returns the index of the value in the range, whether it was added or already existed. \ingroup algos
 * \tparam T The type of the value to be appended.
 * \param r The random access range to be modified.
 * \param v The value to be appended.
 * \return The index of the value in the range.
 */
template <typename T>
auto appendUnique(std::ranges::random_access_range auto& r, T&& v) {
  static_assert(std::is_same_v<typename decltype(begin(r))::value_type, std::remove_reference_t<decltype(v)>>);
  const auto it = find(begin(r), end(r), v);
  size_t index  = std::distance(begin(r), it);
  if (it == end(r))
    r.push_back(std::forward<T>(v));

  return index;
}

/**
 * \brief Prints the contents of a container to the specified output stream.
 *\details
 * This function prints the contents of the given container to the specified output stream.
 * \ingroup algos
 * \tparam C The type of the container to be printed.
 * \param c The container whose contents will be printed.
 * \param os The output stream where the contents will be printed. Default is std::cout.
 */
template <class C>
void printContent(C&& c, std::ostream& os = std::cout) {
  std::ranges::for_each(c, [&os](auto&& var) { os << var << '\n'; });
}

/**
 * \brief Transforms a value range to a pointer range.
 *
 * \ingroup algos
 * \tparam C The type of the container containing values.
 * \param cont The container whose values will be transformed to pointers.
 * \return A subrange containing pointers to the values.
 */
template <class C>
auto transformValueRangeToPointerRange(C& cont) {
  auto transformValueToPointer = [](auto&& obj) { return &obj; };
  return (std::ranges::subrange(cont.begin(), cont.end()) | std::views::transform(transformValueToPointer));
}

/**
 * \brief Transforms a pointer range to a reference range.
 *\details
 * This function transforms a range of pointers to a range of references.
 * \ingroup algos
 * \tparam C The type of the container containing pointers.
 * \param cont The container whose pointers will be transformed to references.
 * \return A subrange containing references to the pointed objects.
 */
template <class C>
auto transformPointerRangeToReferenceRange(C& cont) {
  auto transformValueToPointer = [](auto&& obj) -> auto& { return *obj; };
  return (std::ranges::subrange(cont.begin(), cont.end()) | std::views::transform(transformValueToPointer));
}

#ifndef DOXYGEN
// Forward declare functions
template <typename... Types>
auto makeNestedTupleFlat(std::tuple<Types...> tup);
#endif

namespace Impl {
  template <class Tuple, std::size_t... I>
  constexpr auto makeTupleSubsetImpl(Tuple&& t, std::index_sequence<I...>) {
    return std::make_tuple(std::get<I>(std::forward<Tuple>(t))...);
  }

  template <class Tuple, std::size_t... I>
  constexpr auto makeTupleFromTupleIndicesImpl(Tuple&& t, std::index_sequence<I...>) {
    return std::make_tuple(std::get<I>(std::forward<Tuple>(t))...);
  }

  template <typename T, typename... Ts>
  struct uniqueImpl : std::type_identity<T>
  {
  };

  template <typename... Ts, typename U, typename... Us>
  struct uniqueImpl<std::tuple<Ts...>, U, Us...>
      : std::conditional_t<(std::is_same_v<U, Ts> || ...), uniqueImpl<std::tuple<Ts...>, Us...>,
                           uniqueImpl<std::tuple<Ts..., U>, Us...>>
  {
  };

  template <typename... Ts>
  using unique_tupleImpl = typename uniqueImpl<std::tuple<>, Ts...>::type;

  template <typename T, typename... Types>
  auto makeNestedTupleFlatImpl() {
    constexpr bool isTuple = traits::isSpecialization<std::tuple, T>::value;
    if constexpr (sizeof...(Types) > 0) {
      if constexpr (isTuple)
        return std::tuple_cat(makeNestedTupleFlat(T()), makeNestedTupleFlatImpl<Types...>());
      else
        return std::tuple_cat(std::make_tuple(T()), makeNestedTupleFlatImpl<Types...>());
    } else {
      if constexpr (isTuple)
        return makeNestedTupleFlat(T());
      else
        return std::make_tuple(T());
    }
  }

  template <typename T, typename... Types>
  auto makeNestedTupleFlatAndStoreReferencesImpl(const std::tuple<T, Types...>& tup) {
    constexpr bool isTuple = traits::isSpecialization<std::tuple, std::remove_cvref_t<T>>::value;
    if constexpr (sizeof...(Types) > 0) {
      if constexpr (isTuple)
        return std::tuple_cat(
            makeNestedTupleFlatAndStoreReferencesImpl(std::get<0>(tup)),
            std::apply(
                [](const T&, const Types&... args) {
                  return makeNestedTupleFlatAndStoreReferencesImpl(std::make_tuple(std::cref(args)...));
                },
                tup));
      else
        return std::tuple_cat(
            std::make_tuple(std::cref(std::get<0>(tup))),
            std::apply(
                [](const T&, const Types&... args) {
                  return makeNestedTupleFlatAndStoreReferencesImpl(std::make_tuple(std::cref(args)...));
                },
                tup));
    } else {
      if constexpr (isTuple)
        return makeNestedTupleFlatAndStoreReferencesImpl(std::get<0>(tup));
      else
        return std::make_tuple(std::cref(std::get<0>(tup)));
    }
  }

  template <typename T, typename... Types>
  auto makeNestedTupleFlatAndStoreReferencesNonConstImpl(const std::tuple<T, Types...>& tupconst) {
    auto& tup              = const_cast<std::tuple<T, Types...>&>(tupconst);
    constexpr bool isTuple = traits::isSpecialization<std::tuple, std::remove_cvref_t<T>>::value;
    if constexpr (sizeof...(Types) > 0) {
      if constexpr (isTuple)
        return std::tuple_cat(
            makeNestedTupleFlatAndStoreReferencesNonConstImpl(std::get<0>(tup)),
            std::apply(
                [](T&, Types&... args) {
                  return makeNestedTupleFlatAndStoreReferencesNonConstImpl(std::make_tuple(std::ref(args)...));
                },
                tup));
      else
        return std::tuple_cat(
            std::make_tuple(std::ref(std::get<0>(tup))),
            std::apply(
                [](T&, Types&... args) {
                  return makeNestedTupleFlatAndStoreReferencesNonConstImpl(std::make_tuple(std::ref(args)...));
                },
                tup));
    } else {
      if constexpr (isTuple)
        return makeNestedTupleFlatAndStoreReferencesNonConstImpl(std::get<0>(tup));
      else
        return std::make_tuple(std::ref(std::get<0>(tup)));
    }
  }

} // namespace Impl

/**
 * \brief Finds the index of the first element in the tuple satisfying a predicate.
 *
 * \tparam Tuple Type of the input tuple.
 * \tparam Predicate Type of the predicate function.
 * \param tuple Input tuple.
 * \param pred Predicate function to check each element.
 * \return Index of the first element satisfying the predicate. If no element satisfies the
 *         predicate, it returns the size of the tuple.
 *
 * \details This function takes a tuple and a predicate function and finds the index of the first element in the
 * tuple that satisfies the given predicate. It uses Dune::Hybrid::forEach to iterate through the tuple.
 *
 * \ingroup algos
 */
template <typename Tuple, typename Predicate>
constexpr size_t find_if(Tuple&& tuple, Predicate pred) {
  size_t index        = std::tuple_size<std::remove_reference_t<Tuple>>::value;
  size_t currentIndex = 0;
  bool found          = false;

  Dune::Hybrid::forEach(tuple, [&](auto&& value) {
    if (!found && pred(value)) {
      index = currentIndex;
      found = true;
    }
    ++currentIndex;
  });
  return index;
}

/**
 * \brief Checks if none of the elements in the tuple satisfy a given predicate.
 *
 * \tparam Tuple Type of the input tuple.
 * \tparam Predicate Type of the predicate function.
 * \param tuple Input tuple.
 * \param pred Predicate function to check each element.
 * \return bool True if none of the elements satisfy the predicate, false otherwise.
 *
 *
 * \ingroup algos
 */
template <typename Tuple, typename Predicate>
bool none_of(Tuple&& tuple, Predicate pred) {
  return find_if(tuple, pred) == std::tuple_size<std::decay_t<Tuple>>::value;
}

/**
 * \brief Checks if any of the elements in the tuple satisfy a given predicate.
 *
 * \tparam Tuple Type of the input tuple.
 * \tparam Predicate Type of the predicate function.
 * \param tuple Input tuple.
 * \param pred Predicate function to check each element.
 * \return bool True if any of the elements satisfy the predicate, false otherwise.
 *
 *
 * \ingroup algos
 */
template <typename Tuple, typename Predicate>
bool any_of(Tuple&& tuple, Predicate pred) {
  return !none_of(tuple, pred);
}

/**
 * \brief Filters the elements of a tuple based on a given predicate.
 *
 * \tparam Tuple Type of the input tuple.
 * \tparam Predicate Type of the predicate function.
 * \param tuple Input tuple.
 * \param pred Predicate function to filter the elements.
 * \return auto Tuple containing elements that satisfy the predicate.
 *
 * \details This function applies the given predicate to each element of the tuple. It constructs a new tuple
 * containing only those elements for which the predicate returns true. The resulting tuple is returned.
 *
 * \ingroup algos
 */
template <typename Tuple, typename Predicate>
auto filter(Tuple&& tuple, Predicate pred) {
  return std::apply(
      [&pred](auto... ts) {
        return std::tuple_cat(std::conditional_t<pred(ts), std::tuple<decltype(ts)>, std::tuple<>>{}...);
      },
      tuple);
}

/**
 * \brief Creates a tuple with unique types from the given tuple.
 *
 * \tparam Types Variadic template parameters representing types.
 * \param tuple Input tuple.
 * \return auto Tuple with unique types.
 *
 * \details This function takes a tuple and returns a new tuple containing only unique types from the input tuple.
 *
 * \ingroup algos
 */
template <typename... Types>
constexpr auto unique([[maybe_unused]] std::tuple<Types...>&& tuple) {
  return Impl::unique_tupleImpl<Types...>();
}

/**
 * \brief Counts the number of elements in the tuple satisfying the given predicate.
 *
 * \tparam Tuple Type of the tuple.
 * \tparam Predicate Predicate function determining whether an element satisfies the condition.
 * \param tuple Input tuple.
 * \param pred Predicate function.
 * \return constexpr size_t Number of elements satisfying the predicate.
 *
 * \details This function counts the number of elements in the tuple that satisfy the given predicate.
 *
 * \ingroup algos
 */
template <typename Tuple, typename Predicate>
constexpr size_t count_if(Tuple&& tuple, Predicate pred) {
  size_t counter = 0;
  Dune::Hybrid::forEach(tuple, [&](auto&& value) {
    if (pred(value))
      ++counter;
  });
  return counter;
}

/**
 * \brief Finds the index of the first element in the tuple that is a specialization of the given template type.
 *
 * \tparam Type Template type to search for.
 * \tparam Tuple Type of the tuple.
 * \return Index of the first specialization in the tuple
 *
 * \details This function finds the index of the first element in the tuple that is a specialization of the given
 * template type. It returns the size of the tuple if the template type is not found
 *
 * \ingroup algos
 */
template <template <auto...> class Type, typename Tuple>
constexpr int findTypeSpecialization() {
  return find_if(std::remove_cvref_t<Tuple>(),
                 []<typename T>(T&&) { return traits::isSpecializationNonTypes<Type, std::remove_cvref_t<T>>::value; });
}

/**
 * \brief Gets the specialization of the given template type from the tuple.
 *
 * \tparam Type Template type to search for.
 * \tparam Tuple Type of the tuple.
 * \param tuple The tuple containing elements.
 * \return The specialization element
 *
 * \details This function retrieves the specialization of a template type from the tuple.
 *
 * \ingroup algos
 */
template <template <auto...> class Type, typename Tuple>
auto getSpecialization(Tuple&& tuple) {
  constexpr int index = findTypeSpecialization<Type, Tuple>();
  static_assert(index < std::tuple_size_v<std::remove_cvref_t<Tuple>>,
                "The found index has to be smaller than the tuple size");
  return std::get<index>(tuple);
}

/**
 * \brief Checks if a tuple has a specialization of a template type.
 *
 * \tparam Type Template type to check for.
 * \tparam Tuple Type of the tuple.
 * \return true if the tuple has the specialization; otherwise, false.
 *
 * \details This function checks if a tuple has a specialization of a template type.
 * It uses `find_if` to search for the type and returns true if the index is less than the tuple size; otherwise,
 * false.
 *
 * \ingroup algos
 */
template <template <auto...> class Type, typename Tuple>
constexpr bool hasTypeSpecialization() {
  return (find_if(std::remove_cvref_t<Tuple>(), []<typename T>(T&&) {
            return traits::isSpecializationNonTypes<Type, std::remove_cvref_t<T>>::value;
          }) < std::tuple_size_v<std::remove_cvref_t<Tuple>>);
}

/**
 * \brief Counts the occurrences of a specialization of a template type in a tuple.
 *
 * \tparam Type Template type to count occurrences for.
 * \tparam Tuple Type of the tuple.
 * \return The count of occurrences of the specialization.
 *
 * \details This function counts the occurrences of a specialization of a template type in a tuple.
 *
 * \ingroup algos
 */
template <template <auto...> class Type, typename Tuple>
constexpr bool countTypeSpecialization() {
  return count_if(
      Tuple(), []<typename T>(T&&) { return traits::isSpecializationNonTypes<Type, std::remove_cvref_t<T>>::value; });
}

/**
 * \brief Variable template for counting the occurrences of a specialization of a template type in a tuple.
 *
 * \tparam Type Template type to count occurrences for.
 * \tparam Tuple Type of the tuple.
 *
 * \details This variable template provides a compile-time constant for the count of occurrences
 * of a specialization of a template type in a tuple.
 *
 * \ingroup algos
 */
template <template <auto...> class Type, typename Tuple>
static constexpr bool countTypeSpecialization_v = countTypeSpecialization<Type, Tuple>();

/**
 * \brief Creates a subset tuple with the first N elements from the given tuple.
 *
 * \tparam N Number of elements in the subset.
 * \tparam Tuple Type of the original tuple.
 * \param t The original tuple.
 * \return A new tuple containing the first N elements of the original tuple.
 *
 * \details This function creates a subset tuple with the first N elements from the given tuple.
 *
 * \ingroup algos
 */
template <int N, class Tuple>
constexpr auto makeTupleSubset(Tuple&& t) {
  static_assert(N < std::tuple_size_v<std::remove_reference_t<Tuple>>,
                "The requested size needs to be smaller than the size of the tuple.");

  return Impl::makeTupleSubsetImpl(std::forward<Tuple>(t), std::make_index_sequence<N>{});
}

/**
 * \brief Creates a new tuple using indices from the original tuple.
 *
 * \tparam Tuple Type of the original tuple.
 * \tparam I Indices to include in the new tuple.
 * \param t The original tuple.
 * \return A new tuple containing elements from the original tuple based on the specified indices.
 *
 * \details This function creates a new tuple using indices from the original tuple.
 * It uses `makeTupleFromTupleIndicesImpl` from the `Impl` namespace to implement the tuple creation.
 *
 * \ingroup algos
 */
template <class Tuple, std::size_t... I>
constexpr auto makeTupleFromTupleIndices(Tuple&& t) {
  return Impl::makeTupleFromTupleIndicesImpl(std::forward<Tuple>(t), std::index_sequence<I...>{});
}

/**
 * \brief Creates a flattened nested tuple.
 *
 * \tparam Types Types contained in the original tuple.
 * \return A new flattened nested tuple.
 *
 */
template <typename... Types>
auto makeNestedTupleFlat(std::tuple<Types...>) {
  return decltype(Impl::makeNestedTupleFlatImpl<Types...>())();
}

/**
 * \brief Creates a flattened nested tuple and stores references.
 *
 * \tparam Tuple Type of the original tuple.
 * \param tup The original tuple.
 * \return A new tuple with stored references.
 *
 * \details This function creates a flattened nested tuple and stores references.
 */
template <typename Tuple>
auto makeNestedTupleFlatAndStoreReferences(Tuple&& tup) {
  if constexpr (std::tuple_size_v<std::remove_cvref_t<Tuple>> == 0)
    return tup;
  else if constexpr (!std::is_const_v<std::remove_reference_t<Tuple>>)
    return Impl::makeNestedTupleFlatAndStoreReferencesNonConstImpl(std::forward<Tuple>(tup));
  else
    return Impl::makeNestedTupleFlatAndStoreReferencesImpl(std::forward<Tuple>(tup));
}

/**
 * \brief Returns a reference or std::nullopt if the object is a nullptr.
 * \ingroup tr
 * \tparam T Type of the pointer.
 * \param v Pointer value.
 * \return Reference or std::nullopt.
 */
template <typename T>
requires traits::Pointer<T>
auto& returnReferenceOrNulloptIfObjectIsNullPtr(T v) {
  if constexpr (!std::is_same_v<T, std::nullptr_t>)
    return *v;
  else
    return std::nullopt;
}

} // namespace Ikarus::utils
