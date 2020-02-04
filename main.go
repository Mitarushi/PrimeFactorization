package main

import (
	"fmt"
	"math"
	"math/big"
	"sync"
	"time"
)

type powLog struct {
	log float64
	a   int
	pow int
}

type smooth struct {
	n      int
	factor map[int]int
}

func addSet(a *map[int]*big.Int, n *big.Int) {
	const mod = (1 << 60) - 1
	key := int(new(big.Int).And(n, big.NewInt(mod)).Int64())

	for {
		m, ok := (*a)[key]
		if !ok {
			(*a)[key] = n
			return
		}
		if m.Cmp(n) == 0 {
			return
		}
		key = (key + 1) & mod
	}
}

func addFactor(factor *map[int]int, log float64, cutPows *[7]powLog) {
	if log < 0.1 {
		return
	}
	abs := math.Inf(0)
	minAbs := powLog{0, 0, 0}
	for _, i := range cutPows {
		if math.Abs(log-i.log) < abs {
			abs = math.Abs(log - i.log)
			minAbs = i
		}
	}

	a, ok := (*factor)[minAbs.a]
	if ok {
		(*factor)[minAbs.a] = a + minAbs.pow
	} else {
		(*factor)[minAbs.a] = minAbs.pow
	}
}

func sieve(n, sqrtN *big.Int, x, k int,
	primes *[]int, primeSqrts *[][]int, logPrimes *[]float64,
	smooths chan smooth, wg *sync.WaitGroup) {

	const cutMinLog = 8
	cutMin := big.NewInt(10)
	cutPows := [...]powLog{{0.6931471805599453, 2, 1},
		{1.0986122886681098, 3, 1},
		{1.3862943611198906, 2, 2},
		{1.6094379124341003, 5, 1},
		{1.9459101490553132, 7, 1},
		{2.0794415416798357, 2, 3},
		{2.1972245773362196, 3, 2}}

	table := make([]float64, k)
	factor := make([]map[int]int, k)

	a := new(big.Int).Add(sqrtN, big.NewInt(int64(x)))
	aTemp := new(big.Int).Mul(a, big.NewInt(2))
	aTemp.Add(aTemp, big.NewInt(1))
	b := new(big.Int).Mul(a, a)
	b.Sub(b, n)

	for i := 0; i < k; i++ {
		table[i], _ = new(big.Float).SetInt(b).Float64()
		table[i] = math.Log(table[i])
		b.Add(b, aTemp)
		aTemp.Add(aTemp, big.NewInt(2))

		factor[i] = map[int]int{}
	}

	for i := 0; i < len(*primes); i++ {
		p := big.NewInt(int64((*primes)[i]))
		sqrts := make(map[int]*big.Int)
		for _, v := range (*primeSqrts)[i] {
			addSet(&sqrts, big.NewInt(int64(v)))
		}
		for {
			flag := true
			if p.Cmp(cutMin) == 1 {
				for _, sqrt := range sqrts {
					j := new(big.Int).Sub(sqrt, new(big.Int).Mod(a, p))
					if j.Cmp(big.NewInt(0)) == -1 {
						j.Add(j, p)
					}

					for ; j.Cmp(big.NewInt(int64(k))) == -1; j.Add(j, p) {
						table[j.Int64()] -= (*logPrimes)[i]
						factor[j.Int64()][(*primes)[i]] += 1
						flag = false
					}
				}
			} else {
				flag = false
			}
			if flag {
				break
			}

			newSqrts := make(map[int]*big.Int)
			flag = true
			for _, sqrt := range sqrts {
				temp := new(big.Int)
				temp.Div(temp.Sub(n, temp.Exp(sqrt, big.NewInt(2), nil)), p)
				if (*primes)[i] == 2 {
					if temp.Int64() == 0 {
						addSet(&newSqrts, sqrt)
						addSet(&newSqrts, new(big.Int).Add(p, sqrt))
						addSet(&newSqrts, new(big.Int).Sub(new(big.Int).Mul(p, big.NewInt(2)), sqrt))
						addSet(&newSqrts, new(big.Int).Sub(p, sqrt))
					}
				} else {
					y := new(big.Int).Mul(temp, new(big.Int).ModInverse(new(big.Int).Mul(big.NewInt(2), sqrt), big.NewInt(int64((*primes)[i]))))
					y.Mod(y, big.NewInt(int64((*primes)[i])))
					y.Mul(y, p)
					addSet(&newSqrts, new(big.Int).Add(sqrt, y))
					addSet(&newSqrts, new(big.Int).Sub(p.Mul(p, big.NewInt(int64((*primes)[i]))), new(big.Int).Add(sqrt, y)))
					flag = false
				}
			}
			if flag {
				p.Mul(p, big.NewInt(int64((*primes)[i])))
			}
			sqrts = newSqrts
		}
	}
	count := 0
	for i := 0; i < k; i++ {
		if table[i] > cutMinLog {
			continue
		}
		addFactor(&factor[i], table[i], &cutPows)
		fmt.Printf("%v %v %v\n", i, table[i], factor[i])

		smooths <- smooth{x + i, factor[i]}
		count += 1
	}
	wg.Done()
}

func QS(n *big.Int, b, k int, processNumber, chanStackint int) []int {
	if chanStackint > k {
		chanStackint = k
	}

	if b == -1 {
		nFloat, _ := new(big.Float).SetInt(n).Float64()
		bFloat := math.Exp(math.Sqrt(math.Log(nFloat)*math.Log(math.Log(nFloat))) / 2.0)
		b = int(bFloat) + 1
	}
	fmt.Printf("B is set to %d\n", b)

	start := time.Now()
	rawPrimes := eratosthenes(b)
	fmt.Printf("There are %d primes lower than %d. It took %fs to calculate.\n", len(rawPrimes), b, time.Now().Sub(start).Seconds())

	start = time.Now()
	primes := append(make([]int, 0, len(rawPrimes)/2), 2)
	primeSqrts := append(make([][]int, 0, len(rawPrimes)/2), []int{1})
	logPrimes := append(make([]float64, 0, len(rawPrimes)/2), math.Log(2.0))
	for _, i := range rawPrimes {
		if legendre(n, big.NewInt(int64(i))) == 1 {
			primes = append(primes, i)
			sqrt := int(new(big.Int).ModSqrt(n, big.NewInt(int64(i))).Int64())
			primeSqrts = append(primeSqrts, []int{sqrt, i - sqrt})
			logPrimes = append(logPrimes, math.Log(float64(i)))
		}
	}
	fmt.Printf("There are %d primes in use. It took %fs to sort.\n", len(primes), time.Now().Sub(start).Seconds())

	sqrtN := new(big.Int).Sqrt(n)
	sqrtN.Add(sqrtN, big.NewInt(1))

	smooths := make([]smooth, 0, len(primes)+1)
	wg := sync.WaitGroup{}

	start = time.Now()
	startX := 0
	count := 0
	for {
		smoothsChan := make(chan smooth, processNumber*chanStackint)
		for i := 0; i < processNumber; i++ {
			wg.Add(1)
			go sieve(n, sqrtN, startX, k, &primes, &primeSqrts, &logPrimes, smoothsChan, &wg)
			startX += k
		}
		wg.Wait()
		count += 1
		fmt.Printf("%dth sieve has finished. There are %d B-smooth numbers.\n", count, len(smooths))
		close(smoothsChan)
		for {
			i, ok := <-smoothsChan
			if !ok {
				break
			} else {
				smooths = append(smooths, i)
			}
		}
		if len(smooths) > len(primes) {
			break
		}
	}

	fmt.Printf("There are %d B-smooth numbers. It took %fs to sieve.\n", len(smooths), time.Now().Sub(start).Seconds())

	return nil
}

func main() {
	//n, _ := new(big.Int).SetString("340282366920938463463374607431768211457", 10)
	n, _ := new(big.Int).SetString("1522605027922533360535618378132637429718068114961380688657908494580122963258952897654000350692006139", 10)
	//n, _ := new(big.Int).SetString("147573952589676412927", 10)
	QS(n, 100000000, 1000000000, 4, 100)

}
