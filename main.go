package main

import (
	"fmt"
	"math"
	"math/big"
	"sort"
	"sync"
	"time"
)

type smooth struct {
	n      int
	factor []int
}

func addFactor(factor *[]int, x *big.Int) {
	for _, i := range *factor {
		x.Div(x, big.NewInt(int64(i)))
	}
	xInt := int(x.Int64())
	i := 2
	for {
		if xInt == 1 {
			break
		}
		if xInt%i == 0 {
			xInt /= i
			*factor = append(*factor, i)
		} else {
			i += 1
		}
	}
	sort.Ints(*factor)
}

func calcCutMin(b int, primes *[]int) int {
	pows := make([]int, len(*primes))
	for i := 0; i < len(*primes); i++ {
		pows[i] = 1
	}
	allMulti := 1
	cutMin := 0
	for {
		for i := 0; i < len(*primes); i++ {
			if pows[i]*(*primes)[i] <= cutMin {
				pows[i] *= (*primes)[i]
				allMulti *= (*primes)[i]
			}
		}
		if allMulti > b {
			return cutMin - 1
		}
		cutMin += 1
	}
}

func sieve(n, sqrtN *big.Int, x, k int,
	primes *[]int, primeSqrts *[][]int, logPrimes *[]float32,
	smooths chan smooth, wg *sync.WaitGroup, cutMinInt int, cutMinLog float32) {

	cutMin := big.NewInt(int64(cutMinInt))
	/*cutPows := [...]powLog{{0.6931471805599453, 2, 1},
	{1.0986122886681098, 3, 1},
	{1.3862943611198906, 2, 2},
	{1.6094379124341003, 5, 1},
	{1.9459101490553132, 7, 1},
	{2.0794415416798357, 2, 3},
	{2.1972245773362196, 3, 2}}
	*/
	table := make([]float32, k)
	factor := make([][]int, k)

	a := new(big.Int).Add(sqrtN, big.NewInt(int64(x)))
	aTemp := new(big.Int).Mul(a, big.NewInt(2))
	aTemp.Add(aTemp, big.NewInt(1))
	b := new(big.Int).Mul(a, a)
	b.Sub(b, n)

	for i := 0; i < k; i++ {
		table[i], _ = new(big.Float).SetInt(b).Float32()
		table[i] = float32(math.Log(float64(table[i])))
		b.Add(b, aTemp)
		aTemp.Add(aTemp, big.NewInt(2))

		//factor[i] = make([]int, 0, 3)
	}

	for i := 0; i < len(*primes); i++ {
		p := big.NewInt(int64((*primes)[i]))
		sqrts := make([]*big.Int, 0, 4)
		for _, v := range (*primeSqrts)[i] {
			sqrts = append(sqrts, big.NewInt(int64(v)))
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
						factor[j.Int64()] = append(factor[j.Int64()], (*primes)[i])
						flag = false
					}
				}
			} else {
				flag = false
			}
			if flag {
				break
			}

			newSqrts := make([]*big.Int, 0, 4)
			flag = true
			temp := new(big.Int)
			temp.Div(temp.Sub(n, temp.Exp(sqrts[0], big.NewInt(2), nil)), p)
			if (*primes)[i] == 2 {
				if temp.Int64() == 0 {
					if p.Cmp(big.NewInt(2)) == 0 {
						newSqrts = append(newSqrts, big.NewInt(1))
						newSqrts = append(newSqrts, big.NewInt(3))
					} else {
						newSqrts = append(newSqrts, sqrts[0])
						newSqrts = append(newSqrts, new(big.Int).Add(p, sqrts[0]))
						newSqrts = append(newSqrts, new(big.Int).Sub(new(big.Int).Mul(p, big.NewInt(2)), sqrts[0]))
						newSqrts = append(newSqrts, new(big.Int).Sub(p, sqrts[0]))
					}
				}
			} else {
				y := new(big.Int).Mul(temp, new(big.Int).ModInverse(new(big.Int).Mul(big.NewInt(2), sqrts[0]), big.NewInt(int64((*primes)[i]))))
				y.Mod(y, big.NewInt(int64((*primes)[i])))
				y.Mul(y, p)
				newSqrts = append(newSqrts, new(big.Int).Add(sqrts[0], y))
				newSqrts = append(newSqrts, new(big.Int).Sub(p.Mul(p, big.NewInt(int64((*primes)[i]))), new(big.Int).Add(sqrts[0], y)))
				flag = false
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
		addFactor(&factor[i], new(big.Int).Sub(new(big.Int).Exp(new(big.Int).Add(sqrtN, big.NewInt(int64(x+i))), big.NewInt(2), nil), n))
		fmt.Printf("%v %v %v\n", x+i, table[i], factor[i])

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
	logPrimes := append(make([]float32, 0, len(rawPrimes)/2), float32(math.Log(2.0)))
	for _, i := range rawPrimes {
		if legendre(n, big.NewInt(int64(i))) == 1 {
			primes = append(primes, i)
			sqrt := int(new(big.Int).ModSqrt(n, big.NewInt(int64(i))).Int64())
			primeSqrts = append(primeSqrts, []int{sqrt, i - sqrt})
			logPrimes = append(logPrimes, float32(math.Log(float64(i))))
		}
	}
	fmt.Printf("There are %d primes in use. It took %fs to sort.\n", len(primes), time.Now().Sub(start).Seconds())

	sqrtN := new(big.Int).Sqrt(n)
	sqrtN.Add(sqrtN, big.NewInt(1))

	cutMinInt := calcCutMin(b, &primes)
	fmt.Printf("Factors lower than %d will be cut.\n", cutMinInt)

	smooths := make([]smooth, 0, len(primes)+1)
	wg := sync.WaitGroup{}

	start = time.Now()
	startX := 0
	count := 0
	for {
		smoothsChan := make(chan smooth, processNumber*chanStackint)
		for i := 0; i < processNumber; i++ {
			wg.Add(1)
			go sieve(n, sqrtN, startX, k, &primes, &primeSqrts, &logPrimes, smoothsChan, &wg, cutMinInt, float32(math.Log(float64(b))))
			startX += k
		}
		wg.Wait()
		count += 1
		close(smoothsChan)
		for {
			i, ok := <-smoothsChan
			if !ok {
				break
			} else {
				smooths = append(smooths, i)
			}
		}
		fmt.Printf("%dth sieve has finished. There are %d B-smooth numbers.\n", count, len(smooths))
		if len(smooths) > len(primes) {
			break
		}
	}

	fmt.Printf("There are %d B-smooth numbers. It took %fs to sieve.\n", len(smooths), time.Now().Sub(start).Seconds())

	return nil
}

func main() {
	n, _ := new(big.Int).SetString("340282366920938463463374607431768211457", 10)
	//n, _ := new(big.Int).SetString("1522605027922533360535618378132637429718068114961380688657908494580122963258952897654000350692006139", 10)
	//n, _ := new(big.Int).SetString("147573952589676412927", 10)
	QS(n, 1000000, 1000000, 4, 1000000)

}
