# Praktika
Научная практика по численным методам на C++ и OpenGl выполненная совместно с одногруппником.

Моделирование влияния давления, приложенного к поверхности объектов с различными свойствами.

Вся OpenGl часть взята [отсюда](http://www.songho.ca/opengl/gl_mvc.html).
* Управление камерой wasd+мышь.
* Nx, Ny, Nz - разбиение на тетраэдры.
* cNx, cNy, cNz - разбиение на кубы.
* A, B, C - размер куба по x, y, z.
* Переключаться между кубиками или узлами- стрелки на клавиатуре.
* Пробел - следующий слой.
* Backspace - предыдущий слой.
* Enter - задать свойство E и Mu кубу или узлу.
* Points - привязки узлов и начальные условия.
* Load - загрузить сохраненный объект.
* Save - сохранить.
* P - приложенное давление.
* NextIteration - следующая итерация.

# Test1.txt
Объект на ~4000 тетраэдров.
# Test2.txt
Объект на ~500 тетраэдров.

