## Treść projektu
Przedstawiony algorytm realizuje następujące zagadnienia:
- wczytanie instancji
- usunięcie z wczytanych sekwencji nukleotydowych pozycji o wiarygodności poniżej pewnego założonego progu (próg ten jest ustawiany przez użytkownika)
- utworzenie grafu z wierzchołkami odpowiadającymi wszystkim kilkuliterowym podciągom sekwencji po powyższej operacji (długość podciągów z zakresu od 4 do 9 liter jest drugim parametrem ustawianym przez użytkownika)
- połączenie wierzchołków nieskierowanymi krawędziami, jeśli odpowiadają one takim samym podciągom występującym w różnych sekwencjach, a różnica w pozycjach podciągów wewnątrz sekwencji nie jest większa niż dziesięciokrotność długości podciągu
- wyszukanie w grafie w sposób heurystyczny kliki lub struktury zbliżonej do kliki, w której każda sekwencja wejściowa będzie reprezentowana dokładnie jednym wierzchołkiem
- wypisanie rezultatu na wyjściu w postaci: nr sekwencji wejściowej, nr pozycji w tej sekwencji, dla każdego podciągu wchodzącego w skład kliki (struktury zbliżonej do kliki) oraz sekwencję nukleotydową podciągu
### Format wejściowy
Każda z instancji składa się z dwóch plików, jeden typu fasta zawierający pięć losowo wybranych sekwencji nukleotydowych wraz z ich identyfikatorami, drugi typu qual zawierający oceny wiarygodności odpowiadające wybranym sekwencjom, także z towarzyszącymi im identyfikatorami.
Należy zapewnić w nich obecność motywów (po jednym na instancję) o długości kilkunastu nukleotydów. W obrębie jednej instancji sekwencje powinny zawierać po jednym wystąpieniu motywu, natomiast długość ciągu liczbowego z ocenami wiarygodności musi być ostatecznie taka sama, jak długość odpowiedniej sekwencji nukleotydowej.
### Informacje o implementacji
W programie zostały utworzone dwie struktury Data oraz Graph. Pierwsza z nich przechowuje informacje dotyczące sekwencji wejściowych, tj. identyfikator sekwencji,
pozycje oraz ciąg nukleotydów. Druga natomiast reprezentuje graf nieskierowany. Składa się z listy sąsiedztwa i wektora przechowującego ilość sąsiadów dla każdego wierzchołka.
Metoda void addVertex w niej zawarta umożliwia dodawanie wierzchołków do grafu i utworzenie odpowiednich połączeń.

W pierwszej kolejności wczytywane są dane z obu plików fasta oraz qual za pomocą funkcji
INPUT, wyodrębniane są poszczególne informacje (identyfikator, ciąg nukleotydów,
wiarygodność). Dodatkowo następuje weryfikacja wprowadzonych danych pod kątem
zgodności ilości nukleotydów do ilości ocen wiarygodności w każdej sekwencji. Zgodnie z
podaną przez użytkownika minimalną wartością oceny jakości nukleotydu, wszelkie odczyty
mające mniejszą wiarygodnością są usuwane za pomocą funkcji Remove.

Użytkownik na wstępie jest proszony o podanie informacji dotyczącej wielkości
oligonukleotydów (podciągów), która zostanie użyta w celu utworzenia wierzchołków. Na
podstawie danych o długości wszystkich sekwencji po redukcji oraz wielkości podciągów
ustalana jest liczba wierzchołków w grafie. Odpowiedzialna za inicjalizację grafu jest funkcja
Graph_maker. 
### Wyszukiwanie maksymalnej kliki
Przedstawione funkcje do wyszukiwania kliki w grafie odpowiadającej
występowaniu motywu zostały zaimplementowane w oparciu o treści pochodzące z artykułu:
“Fast Algorithms for the Maximum Clique Problem on Massive Sparse Graphs”, autorstwa
Bharath Pattabiraman, Md. Mostofa Ali Patwary, Assefaw H. Gebremedhin, Wei-keng Liao i
Alok Choudhary.

Poszukiwania kliki w tym podejściu heurystycznym opiera się na procedurze rekurencyjnej,
wykorzystującej funkcje CliqueHeu oraz MaxCliqeHeu. Argumentem, który trafia do
procedury jest obiekt struktury Graph. To podejście zostało oparte na rozpatrywaniu tylko wierzchołków o najwyższym stopniu,
konsekwencją tego jest przyspieszenie obliczeń z jednoczesnym ubytkiem na jakości
rozwiązań.
### Złożoność obliczeniowa
Najważniejsze z punktu widzenia teorii złożoności obliczeniowej w przedstawionym
algorytmie są funkcje poszukujące maksymalnej kliki. Algorytm iteruje po n wierzchołkach,
za każdym razem wywołując funkcję rekurencyjną CliqueHeu, która wykonuje się aż do
momentu, gdy podany zbiór będzie pusty. Jest to oczywiście zależne od maksymalnego
stopnia wierzchołka ( ∆ ). Funkcja rekurencyjna wyznacza również sąsiadów wierzchołka, co
wpływa na złożoność całego programu. Jest to zatem algorytm heurystyczny o złożoności wielomianowej.



