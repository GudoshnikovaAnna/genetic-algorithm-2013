using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GeneticApproximation
{
    class Program
    {
        #region declaration
        /*Входные параметры*/
    
        /*Количество хромосом в геноме - тобишь коэффициенты в полиноме*/
        private static int cromosomCount = 5;
        private static double[] genom;  //массив коэффициентов
        private static List<double[]> GenomOldList; // массив геномов
        private static List<double[]> GenomCrossList; // массив геномов после скрещивания
        private static List<double[]> GenomMutantList; // массив геномов после мутации
        private static double[] PolynFunctionList;
        private static double[,] ApproxFuncCoordinates;
        private static int PairsCoordinateCount; // сколько точек подаются на вход
        private static Random rnd = new Random();
 
        private static List<double[]> fitnessResult; //Список из массивов, где первые 5 элементов - геном, последний - значение фитнес-функции
        /*Максимальное количество поколений (итераций)*/
        int maxGeneration = 100;

        /*Максимальный размер популяции - геномов в популяции*/
        private static int populationMaxSize = 20;
        /*Вероятность мутации*/
        double MutationPosibility = 0.1;
        /*Вероятность скрещивания*/
        double CrossoverPosibility = 0.9;
        /*Количество особей, отбираемых в каждом поколении*/
        static int SelectCount = 10;
        #endregion 

        /*Функция для получения матрицы координат апроксимируемой функции*/
        private static double [,] GetApproxFuncCoordinates()
        {
            /*Просим пользователя подать на вход текстовый документ с координатами апроксимуремой функции*/
            Console.WriteLine("Please, enter file name with coordinates:");
            string fileName = Console.ReadLine();
            StreamReader sr = File.OpenText(@fileName);
            string matrix = sr.ReadToEnd(); // здесь матрица в строковом виде
            /*считаем размер матрицы с координатами, знаю криво, но пока пусть так*/
            PairsCoordinateCount = 0; // PairsCoordinateCount - число пар координат
            
            for (int i = 0; i < matrix.Length; i++)
            {
                if (matrix[i].Equals('\n')) PairsCoordinateCount++;
            }
            PairsCoordinateCount++;
            Console.WriteLine("Количество пар " + PairsCoordinateCount);
            
            matrix = matrix.Replace("\r\n", "");
            /*Заполняем матрицу*/
            double[,] ApproxFuncCoordinates = new double[PairsCoordinateCount, 2];

            double[] list = new double[PairsCoordinateCount * 2];
            int h = 0;
            for (int i = 0; i < matrix.Length; i++)
            {
                /*StringBuilder - изменяющая строка, т.е. почти как строковый буфер*/
                StringBuilder sb = new StringBuilder();
                /*короче, накапливаем в этом буфере число до ',' 
                 * затем переписываем в список, уже преобразованное в double*/
                while (!(matrix[i].Equals(';')))
                {
                    sb.Append(matrix[i++]);
                }
                list[h++] = double.Parse(sb.ToString());// list [x1,y1,x2,y2, ...]
            }
            PolynFunctionList = new double[PairsCoordinateCount];
            /*Заполняем нашу матрицу числами из списка*/
            for (int i = 0; i < PairsCoordinateCount; i++)
                for (int j = 0; j < 2; j++)
                    ApproxFuncCoordinates[i, j] = list[j + i * 2]; 
            return ApproxFuncCoordinates;

        }
        //Начальный геном заполняем рандомно
        private static double[] getGenom_start()
        {
           genom = new double[cromosomCount];
           Console.WriteLine("genom:");
            for (int i = 0; i < genom.Length; i++)
            {
               genom[i] = rnd.Next(-500, 500);
               Console.WriteLine(genom[i]);
            }
            Console.WriteLine();
            return genom;
        }

        //Считаем полином по всем данным х-ам (которые подавались на вход), со сгенеринными коэффициентами
        private static void getPolynFunctionList(double[] genom)
        {             
            for (int j = 0; j < PairsCoordinateCount; j++)
            {
                PolynFunctionList[j] = genom[0] * Math.Pow(ApproxFuncCoordinates[j, 0], 4) + 
                    genom[1] * Math.Pow(ApproxFuncCoordinates[j, 0], 3) +
                    genom[2] * Math.Pow(ApproxFuncCoordinates[j, 0], 2) + 
                    genom[3] * ApproxFuncCoordinates[j, 0] + genom[4];
                Console.WriteLine( "f = "+ PolynFunctionList[j]);
            }     
        }
        
        private static void GenericStartPopulation()
        {
            GenomOldList = new List<double[]>();
            //fitnessResult = new List<double[]>(); //список из массивов=геном+фитфункция
            Dictionary<double, double[]> FitnessGenom = new Dictionary<double, double[]>();
            double[] fitvalue = new double[populationMaxSize]; //значения фитфункции
            double[] genomFitValue = new double[cromosomCount + 1]; //геном + значение фитнес-функции
            for (int i = 0; i < populationMaxSize; i++)
            {
                getGenom_start(); //генерируем популяцию
                GenomOldList.Add(genom);
                getPolynFunctionList(genom); //считаем значения полинома с коэффициентами-хромосомами генома
                FitnessGenom[fitnessFunction()] = genom;
                fitvalue[i] = fitnessFunction();
            }

            Array.Sort(fitvalue); //сортируем значения фитнес-функций по возрастанию
            double[] SelectedFitValue = new double[SelectCount]; //массив с топ-10 значениями фитфункции
            for (int k = 0; k < SelectCount; k++)
            {
                SelectedFitValue[k] = fitvalue[k];
            }

            List<double[]> bestGenom = new List<double[]>();//геномы с лучшими фитнес-функциями
            foreach (var elem in FitnessGenom)
            {
                if (SelectedFitValue.Contains(elem.Key))
                {
                    bestGenom.Add(elem.Value);
                }
            }
        }

        /*Надо посчитать фитнесс функцию - фукнцию приспособленности. 
         * Сначала сосчитаем разницу между истинным значанием функции 
         * (т.е. с тем, что подаем на вход) и с тем, что получили из популяции
         * затем берем квадрат этой разности
         * и складывем
         * Т.о. определили метрику
         */
        private static double fitnessFunction() 
        {
            double sum = 0;
            for (int i = 0; i < PairsCoordinateCount; i++)
            {
                sum = sum + Math.Pow((PolynFunctionList[i] - ApproxFuncCoordinates[i,1]), 2); 
            }
            Console.WriteLine("sum " + sum);
            return sum; 
        }// надо ее уменьшить. 

        // скрещивание от двух родителей
        private static double[] GenomCrossover(double[] parentGenom1, double[] parentGenom2)
        {
            double[] childGenom = new double[cromosomCount];
           
            for (int i = 0; i < cromosomCount; i++)
            {
                childGenom[i] = (parentGenom1[i] + parentGenom2[i]) / 2;
            }
            return childGenom;
        }

        

        // новое поколение геномов, записанные в список GenomCrossList после скрещивания
        private static void NewCrossoverPopulation()
        {
            GenomCrossList = new List<double[]>();
            double[] crossgenom = new double[cromosomCount];
            double[] testArr = new double[5];
            testArr = GenomOldList[1];
            Console.WriteLine("GenomOldList");
            for (int t = 0; t < 5; t++)
            {
                Console.Write(" " + testArr[t]);
            }
            for (int i = 0; i < populationMaxSize-1; i++)
                for (int j = i + 1; j < populationMaxSize; j++)
                {
                    
                    //crossgenom = GenomCrossover(GenomOldList[i],GenomOldList[j]);
                    //GenomCrossList.Add(crossgenom);
                }
            //Console.Out.WriteLine(GenomCrossList.Count());
        }

        // мутация генома
        private static double[] GenomMutation(double[] crossgenom)
        {
            int chrom = rnd.Next(0, crossgenom.Length-1);
            crossgenom[chrom] = rnd.Next(-500, 500);
            return crossgenom;
        }
        //популяция мутировавших геномов
        private static void NewPopulation()
        {
            for (int i = 0; i < 29; i++)
            {
                int k = rnd.Next(0, populationMaxSize-1);

            }

                foreach (double[] crossgenom in GenomCrossList)
                {
                    double[] mutationgenom = GenomMutation(crossgenom);
                    GenomMutantList.Add(mutationgenom);
                }
        }

        private static void Selection()
        {
            
            
        }

        static void Main(string[] args)
        {
            ApproxFuncCoordinates = GetApproxFuncCoordinates();
            GenericStartPopulation();
           //Console.WriteLine("Count old list " + GenomOldList.Count());
           // GenericPopulation_Start();
            //NewCrossoverPopulation();
            Console.WriteLine("Finish");
            Console.Read();
        }
    }
}
