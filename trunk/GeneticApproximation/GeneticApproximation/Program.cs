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
        //private static double[] PolynFunctionList;
        private static double[,] ApproxFuncCoordinates;
        private static int PairsCoordinateCount; // сколько точек подаются на вход
        private static Random rnd = new Random();
        private static Dictionary<double, double[]> FitnessGenom = new Dictionary<double, double[]>();//ключ - значение фитнес-функции, значение - соответствующий геном
        private static double MinFitnessFunction;
        /*Максимальный размер популяции - геномов в популяции*/
        private static int populationMaxSize = 10;
        /*Вероятность скрещивания*/
        private static double CrossoverPosibility = 0.9;
        /*Количество особей, отбираемых в каждом поколении*/
        static int SelectCount = 7;
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
            //Console.WriteLine("Количество пар " + PairsCoordinateCount);
            
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
            //PolynFunctionList = new double[PairsCoordinateCount];
            /*Заполняем нашу матрицу числами из списка*/
            for (int i = 0; i < PairsCoordinateCount; i++)
                for (int j = 0; j < 2; j++)
                    ApproxFuncCoordinates[i, j] = list[j + i * 2]; 
            return ApproxFuncCoordinates;

        }
        
        //Начальный геном заполняем рандомно
        private static double[] getGenom_start()
        {
           double[] genom = new double[cromosomCount];
          // Console.WriteLine("genom:");
            for (int i = 0; i < genom.Length; i++)
            {
               genom[i] = rnd.Next(-500, 500);
              // Console.WriteLine(genom[i]);
            }
            //Console.WriteLine();
            return genom;
        }

        //Считаем полином по всем данным х-ам (которые подавались на вход), со сгенеринными коэффициентами
        private static double[] getPolynFunctionList(double[] genom)
        {
            double[] PolynFunctionList = new double[PairsCoordinateCount];
            for (int j = 0; j < PairsCoordinateCount; j++)
            {
                PolynFunctionList[j] = genom[0] * Math.Pow(ApproxFuncCoordinates[j, 0], 4) + 
                    genom[1] * Math.Pow(ApproxFuncCoordinates[j, 0], 3) +
                    genom[2] * Math.Pow(ApproxFuncCoordinates[j, 0], 2) + 
                    genom[3] * ApproxFuncCoordinates[j, 0] + genom[4];
               // Console.WriteLine( "f = "+ PolynFunctionList[j]);
            }
            return PolynFunctionList;
        }
        
        private static List<double[]> GenerateStartPopulation()
        {
            List<double[]> StartPopulation = new List<double[]>();
            double[] genom;
            for (int i = 0; i < populationMaxSize; i++)
            {
                genom = getGenom_start(); //генерируем геном
                StartPopulation.Add(genom);
            }
            return StartPopulation;
        }

        private static List<double[]> getGenomsForSelection(List<double[]> population)
        {
            FitnessGenom.Clear();
            List<double[]> bestGenom = new List<double[]>();//геномы с лучшими фитнес-функциями            
            List<double> tempFitFuncValue = new List<double>(); //временный список, чтобы сохранить значения фитнес-функций

            foreach (var genom in population)
            {
                //getPolynFunctionList(genom); //посчитали полиномы с коэффициентами-хромосомами для всех входных иксов
                double temp = fitnessFunction(genom); //посчитали фитнес-функцию для текущего генома
                tempFitFuncValue.Add(temp);//добавили значение в список
                FitnessGenom[temp] = genom;//добавили значение в словарь
            }
            double[] fitvalue = new double[tempFitFuncValue.Count]; //значения фитфункции
            /*переносим список в массив, чтобы быстро автоматически отсортировать потом*/
            for (int i = 0; i < tempFitFuncValue.Count; i++)
            {
                fitvalue[i] = tempFitFuncValue[i];
            }
            
            Array.Sort(fitvalue); //сортируем значения фитнес-функций по возрастанию
            double[] SelectedFitValue = new double[SelectCount]; //массив с топ-10 значениями фитфункции
            for (int k = 0; k < SelectCount; k++)
            {
                SelectedFitValue[k] = fitvalue[k];
                //Console.WriteLine("после отсеивания = " + SelectedFitValue[k]);
            }

            MinFitnessFunction = SelectedFitValue[0];

            foreach (var elem in FitnessGenom)
            {
                for (int p = 0; p < SelectedFitValue.Length; p++)
                {
                    if (elem.Key == SelectedFitValue[p])
                    {
                        bestGenom.Add(elem.Value);
                    }
                }
            }
           
            return bestGenom; 
        }

        /*Надо посчитать фитнесс функцию - фукнцию приспособленности. 
         * Сначала сосчитаем разницу между истинным значанием функции 
         * (т.е. с тем, что подаем на вход) и с тем, что получили из популяции
         * затем берем квадрат этой разности
         * и складывем
         * Т.о. определили метрику
         */
        private static double fitnessFunction(double[] genom) 
        {
            double sum = 0;
            double[] PolynFunctionList = getPolynFunctionList(genom); 
            for (int i = 0; i < PairsCoordinateCount; i++)
            {
                sum = sum + Math.Pow((PolynFunctionList[i] - ApproxFuncCoordinates[i,1]), 2); 
            }
            return sum; 
        }// надо ее уменьшить. 

        // скрещивание от двух родителей
        private static double[] GenomCrossover(double[] parentGenom1, double[] parentGenom2)
        {
            double[] childGenom = new double[cromosomCount];
            
            for (int i = 0; i < cromosomCount; i++)
            {
                {
                    childGenom[i] = (parentGenom1[i] + parentGenom2[i]) / 2;
                }
            }
            return childGenom;
        }

        // новое поколение геномов, записанные в список GenomCrossList после скрещивания        
        private static List<double[]> NewCrossoverPopulation(List<double[]> OldPopulation)
        {
            List<double[]> selectPopulation =  getGenomsForSelection(OldPopulation);
            double[] newgenom = new double[cromosomCount];
            List<double[]> NewPopulation = new List<double[]>();
            for (int i = 0; i < selectPopulation.Count - 1; i++)
            {
                for (int j = i + 1; j < selectPopulation.Count; j++)
                {
                    double p = rnd.NextDouble();
                    if (p < CrossoverPosibility)
                    {
                        newgenom = GenomCrossover(selectPopulation[i], selectPopulation[j]);
                        NewPopulation.Add(newgenom);
                    }
                    else
                    {
                        NewPopulation.Add(selectPopulation[i]);
                    }
                }
            }
            return NewPopulation;
        }

        // мутация генома
        private static double[] GenomMutation(double[] crossgenom)
        {
            int chrom = rnd.Next(0, crossgenom.Length-1); //ген (хромосома), который будет мутировать
            crossgenom[chrom] = rnd.Next(-500, 500);
            
            //Console.WriteLine("mutant gen " + crossgenom[chrom]);
            //Console.WriteLine();
            return crossgenom;
        }

        //популяция мутировавших геномов
        private static List<double[]> NewPopulation(List<double[]> OldPopulation)
        {
            double[] Mutantgenom;//геном, который мутирует
            List<double[]> XMen = new List<double[]>();//окончательное новое поколение после скрещивания и мутации
            List<double[]> CrossPop = NewCrossoverPopulation(OldPopulation);
            XMen = CrossPop;

            int MutantNumber = 3;//rnd.Next(1,XMen.Count - 1); // сколько мутирует
            int[] mutantIndex = new int[MutantNumber]; // индексы мутировавших геномов

            for (int i = 0; i < MutantNumber; i++)
                mutantIndex[i] = -1;

            //Console.WriteLine();
            //Console.WriteLine("Осторожно, мутанты отаке");
            for (int i = 0; i < MutantNumber; i++)
            {
                //double p = rnd.NextDouble(); //вероятность мутации
                //if (p > 0.6)
                //{
                    int k = rnd.Next(0, CrossPop.Count - 1); //индекс мутирующего генома
                    if (!mutantIndex.Contains(k)) //проверяем, не мутировал ли уже этот геном, если нет - велкам
                    {
                        Mutantgenom = GenomMutation(CrossPop[k]);                       
                        mutantIndex[i] = k;
                        XMen.RemoveAt(k); //удаляем старый геном
                        XMen.Insert(k, Mutantgenom); //вписываем новый, мутировавший
                        //getPolynFunctionList(Mutantgenom);
                        //Console.WriteLine("фитнесс-функция от мутировавшего генома №" + k +" = " + fitnessFunction(Mutantgenom));
                    }
                //}
            }
            return XMen;
        }

        private static void GeneticAlg()
        {
            List<double[]> TempPopulation = new List<double[]>();
            double[] FinalGenom;
            TempPopulation = NewPopulation(GenerateStartPopulation());
            //Console.WriteLine("MinFitness = " + MinFitnessFunction);
            int count = 0;
            while (count < 100000)
            {
               FitnessGenom.Clear();
               List<double[]> temp =  NewPopulation(TempPopulation);
                /*for (int k = 0; k < temp.Count; k++)
                {
                    if (k >= populationMaxSize)
                        temp.RemoveAt(k);
                }*/
                TempPopulation = temp;
               
                //Console.WriteLine("MinFitness = " + MinFitnessFunction);
                count++;
            }
            FinalGenom = FitnessGenom[MinFitnessFunction];

            Console.WriteLine();
            Console.WriteLine("FinalGenom: ");
            for (int k = 0; k < cromosomCount; k++)
            {
                Console.WriteLine(FinalGenom[k]);
            }
            Console.WriteLine();
            Console.WriteLine("MinFitness = " + MinFitnessFunction);
        }

        static void Main(string[] args)
        {
            ApproxFuncCoordinates = GetApproxFuncCoordinates();
            GeneticAlg();        
            Console.WriteLine("Finish");
            Console.Read();
        }
    }
}