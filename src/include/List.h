/*! \file List.h
 *
 * \brief A library of classes used to manage dynamic lists.
 */

#ifndef INCLUDE_LIST_H_
#define INCLUDE_LIST_H_

/*! \brief A class to represent dynamic lists.
 *
 * This class helps in the creation and management of dynamic lists of
 * objects, whose size is not known in advance. List offers the advantage
 * of saving memory, since only the necessary space will be allocated,
 * but it has the disadvantage of creating an object which is not contiguous
 * in memory and, therefore, is very inefficient for subsequent manipulation.
 *
 * For this reason, the best use of List objects is to collect all the
 * desired members and then, once the element number is known, to convert
 * the List to C array, by calling `List.to_array()`. This function returns
 * a contiguous array of type T[SIZE] that can be used for indexed access.
 */
template<class T> class List {
protected:
  int size; //!< Size of the List.
  struct element {
    T value; //!< Value of the list element.
    element* p_prev; //!< Pointer to the previous element in the list.
  }; //!< List element connector.
  element *current, //!< Pointer to element affected by last operation.
    *first, //!< Pointer to the first element in the List.
    *last; //!< Pointer to the last element in the List.

public:
  /*! \brief List constructor.
   *
   * Use the constructor `List<T>([int length])` to create a new list with a given
   * size. If the required size is not known in advance, it is recommended
   * to create a List with SIZE=1 (this is the default behavior) and then
   * to append the elements dynamically, using `List.append(ELEMENT)` (where
   * ELEMENT needs to be a value of type T, corresponding to the class
   * template specialization). Note that, due to the default behavior, the
   * following calls are equivalent and they both produce an integer List
   * with size equal to 1:
   *
   * \code{.cpp}
   * List<int> a = List<int>(1);
   * List<int> b = List<int>();
   * \endcode
   *
   * \param length: `int` The size of the list to be constructed [OPTIONAL, default=1].
   */
  List(int length = 1) {
    size = length;
    first = new element;
    first->p_prev = NULL;
    current = first;
    element *p_prev = first;
    for (int i = 1; i < size; i++) {
      current = new element;
      current->p_prev = p_prev;
      p_prev = current;
    }
    last = current;
  }

  /*! \brief Destroy a List instance.
   */
  ~List() {
    current = last;
    element *old;
    while (current->p_prev) {
      old = current;
      current = old->p_prev;
      delete old;
    }
    delete current;
  }
  
  /*! \brief Append an element at the end of the List.
   *
   * To dynamically create a list whose size is not known in advance,
   * elements can be appended in an iterative way. Note that element
   * manipulation is much more effective in a C array than in a List
   * object. For this reason, after the List has been created, it is
   * strongly advised to convert it to a C array by calling the function
   * `List.to_array()`.
   *
   * \param value: `T` The value of the element to be appended.
   */
  void append(T value) {
    element *p_prev = last;
    current = new element;
    current->value = value;
    current->p_prev = p_prev;
    last = current;
    size++;
  }

  /*! \brief Get the element at given index.
   *
   * Get the element specified by the index argument. The first element
   * has index 0 and the last one has index [size - 1].
   *
   * \param index: `int` The index of the element to be retrieved. 0 for first.
   * \return value `T` The value of the element at the requested position.
   * \throws ListOutOfBoundsException: Raised if the index is out of bounds.
   */
  T get(int index) {
    if (index < 0 || index > size - 1) {
      throw ListOutOfBoundsException(index, 0, size - 1);
    }
    current = last;
    for (int i = size - 1; i > index; i--) current = current->p_prev;
    return current->value;
  }

  /*! \brief Get the number of elements in the List.
   *
   * Get the number of elements currently stored in the List.
   *
   * \return size `int` The size of the List.
   */
  int length() {
    return size;
  }

  /*! \brief Set an element by index and value.
   *
   * Set the element at the position specified by the index to the value
   * specified by the value argument.
   *
   * \param index: `int` The index of the element to be set. 0 for first.
   * \param value: `int` The value to store in the pointed element.
   * \throws ListOutOfBoundsException: Raised if the index is out of bounds.
   */
  void set(int index, T value) {
    if (index < 0 || index > size - 1) {
      throw ListOutOfBoundsException(index, 0, size - 1);
    }
    current = last;
    for (int i = size - 1; i > index; i--) current = current->p_prev;
    current->value = value;
  }

  /*! \brief Convert the list to a C array.
   *
   * The List object is useful to dynamically manage a set of objects
   * when the number of elements is not known in advance. However, the
   * resulting object is not contiguosly stored in memory. As a result,
   * access to specific elements in the middle of the list is not very
   * effective, because the list needs to be walked every time up to
   * the desired position. In order to avoid this, `List.to_array()` makes
   * a conversion from List to C array, returning a contiguous object,
   * where indexed access can be used.
   *
   * \return array `T*` A C array of type T and size equal to the List size.
   */
  T* to_array() {
    T *array = new T[size];
    current = last;
    for (int i = size - 1; i > -1; i--) {
      array[i] = current->value;
      current = current->p_prev;
    }
    return array;
  }
};

#endif
