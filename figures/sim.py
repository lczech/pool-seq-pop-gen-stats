import msprime
import numpy as np
import pandas as pd


# num_x_mat is shape [num_SNPs][num_pops]
def print_count_matrix(num_0_mat, num_1_mat, outfile):
    assert num_0_mat.shape == num_1_mat.shape
    drop = np.logical_or(np.all(num_0_mat == 0, axis=1),
                         np.all(num_1_mat == 0, axis=1))
    num_0_mat = num_0_mat[~drop, :]
    num_1_mat = num_1_mat[~drop, :]
    with open(outfile, 'w') as ofh:
        for i in range(len(num_0_mat)):
            gt_strs = []
            for p in range(num_0_mat.shape[1]):
                gt_strs.append(
                    '{}:{}:0:0:0:0'.format(num_0_mat[i, p],
                                           num_1_mat[i, p])
                )
            ofh.write(
                '\t'.join(['chr1', str(i+1), 'A'] + gt_strs) + '\n'
            )


def sim_sequencing(pool_freqs, avg_read_depth, seq_error):
    read_depths = np.random.poisson(avg_read_depth, size=pool_freqs.shape)
    true_1_counts = np.random.binomial(read_depths, pool_freqs)
    true_0_counts = read_depths - true_1_counts
    swap_1_to_0 = np.random.binomial(true_1_counts, seq_error)
    swap_0_to_1 = np.random.binomial(true_0_counts, seq_error)

    obs_1_counts = true_1_counts - swap_1_to_0 + swap_0_to_1
    obs_0_counts = true_0_counts - swap_0_to_1 + swap_1_to_0

    assert np.all(obs_1_counts + obs_0_counts == read_depths)

    return obs_0_counts, obs_1_counts


def statistical_model(true_freqs, pool_size, avg_read_depth, seq_error):
    pool_freqs = np.random.binomial(pool_size, true_freqs) / pool_size
    return sim_sequencing(pool_freqs, avg_read_depth, seq_error)


# true freqs is shape [num_SNPs][num_pops]
def theta_pi(true_freqs):
    assert true_freqs.shape[1] == 1
    return np.mean(2 * true_freqs * (1-true_freqs))


def pi_within(true_freqs):
    assert true_freqs.shape[1] == 2

    single_pop_pi = np.mean(2 * true_freqs * (1-true_freqs),
                            axis=0)
    return np.mean(single_pop_pi)


def pi_between(true_freqs):
    assert true_freqs.shape[1] == 2

    p0 = np.mean((1-true_freqs[:, 0])*true_freqs[:, 1])
    p1 = np.mean(true_freqs[:, 0]*(1-true_freqs[:, 1]))
    return p0 + p1


def pi_total(true_freqs):
    assert true_freqs.shape[1] == 2

    avg_freqs = np.mean(true_freqs, axis=1)
    return np.mean(2 * avg_freqs * (1-avg_freqs))


def fst_hudson(true_freqs):
    assert true_freqs.shape[1] == 2

    pi_b = pi_between(true_freqs)
    pi_w = pi_within(true_freqs)
    return (pi_b - pi_w) / pi_b


def fst_nei(true_freqs):
    assert true_freqs.shape[1] == 2

    pi_t = pi_total(true_freqs)
    pi_w = pi_within(true_freqs)
    return (pi_t - pi_w) / pi_t


def msprime_model(pool_size,
                  read_depth,
                  seq_err,
                  theta,
                  seq_len,
                  num_pops,
                  d_time=None,
                  exponential_growth=False):

    samples = {}
    demo = msprime.Demography()
    if exponential_growth:
        for p in range(num_pops):
            demo.add_population(name='pop' + str(p),
                                initial_size=np.exp(5),
                                growth_rate=10.)
            samples['pop' + str(p)] = pool_size + 1000
        for p in range(1, num_pops):
            demo.add_population_split(d_time + p*1e-10,
                                      derived=['pop'+str(p)],
                                      ancestral='pop0')
    else:
        for p in range(num_pops):
            demo.add_population(name='pop' + str(p),
                                initial_size=1.)
            samples['pop' + str(p)] = pool_size + 1000
        for p in range(1, num_pops):
            demo.add_population_split(d_time + p*1e-10,
                                      derived=['pop' + str(p)],
                                      ancestral='pop0')
    print('Simulating Ancestry')
    ts = msprime.sim_ancestry(samples=samples,
                              demography=demo,
                              ploidy=1,
                              recombination_rate=0.0001,
                              sequence_length=int(seq_len))
    print('Adding mutations')
    mts = msprime.sim_mutations(ts,
                                rate=theta,
                                model=msprime.BinaryMutationModel())
    print('Parsing', mts.num_sites, 'variants')
    p1_pool_freqs = []
    p1_true_freqs = []
    if num_pops == 2:
        p2_pool_freqs = []
        p2_true_freqs = []
    for v in mts.variants():
        gts = v.genotypes
        p1_pool_freqs.append(gts[0:pool_size].mean())
        p1_true_freqs.append(gts[0:(pool_size+1000)].mean())
        if num_pops == 2:
            p2_pool_freqs.append(
                gts[(pool_size+1000):(2*pool_size+1000)].mean()
            )
            p2_true_freqs.append(
                gts[(pool_size+1000):(2*pool_size+2000)].mean()
            )
    if num_pops == 1:
        pool_freqs = np.zeros((int(seq_len), 1))
        pool_freqs[0:len(p1_pool_freqs), 0] = p1_pool_freqs
        true_freqs = np.zeros((int(seq_len), 1))
        true_freqs[0:len(p1_true_freqs), 0] = p1_true_freqs
        num_seg = np.sum(pool_freqs > 0)
        return sim_sequencing(pool_freqs,
                              read_depth,
                              seq_error), true_freqs, num_seg
    if num_pops == 2:
        pool_freqs = np.zeros((int(seq_len), 2))
        pool_freqs[0:len(p1_pool_freqs), 0] = p1_pool_freqs
        pool_freqs[0:len(p1_pool_freqs), 1] = p2_pool_freqs
        true_freqs = np.zeros((int(seq_len), 2))
        true_freqs[0:len(p1_true_freqs), 0] = p1_true_freqs
        true_freqs[0:len(p1_true_freqs), 1] = p2_true_freqs
        return sim_sequencing(pool_freqs,
                              read_depth,
                              seq_error), true_freqs


res = pd.DataFrame({'fname': [],
                    'true_pairwise_het': [],
                    'true_num_seg': [],
                    'true_pi_w': [],
                    'true_pi_b': [],
                    'true_pi_t': [],
                    'true_hudson_fst': [],
                    'true_nei_fst': [],
                    'pool_size': [],
                    'read_depth': [],
                    'seq_error': [],
                    'num_sites': []})

for pool_size in [10, 100, 1000]:
    for read_depth in [10, 100, 1000]:
        for seq_error in [0., 1e-5, 1e-4, 1e-3]:
            for rep in range(10):
                print(pool_size, read_depth, seq_error, rep)
                print('stat 1 pop')
                # statisitcal 1 pop
                for alpha in np.linspace(0.1, 0.4, num=5):
                    true_freqs = np.random.beta(alpha, 0.5,
                                                size=int(1e6))
                    true_freqs = true_freqs.reshape((-1, 1))
                    n0, n1 = statistical_model(true_freqs,
                                               pool_size,
                                               read_depth,
                                               seq_error)
                    fname = '{}_{}_{}_{}_stat_1_pop_alpha_{}.txt'.format(
                        pool_size, read_depth, seq_error, rep, alpha
                    )
                    print_count_matrix(n0, n1, fname)
                    res = res.append(
                        {'fname': fname,
                         'true_pairwise_het': theta_pi(true_freqs),
                         'true_num_seg': np.nan,
                         'true_pi_w': np.nan,
                         'true_pi_b': np.nan,
                         'true_pi_t': np.nan,
                         'true_hudson_fst': np.nan,
                         'true_nei_fst': np.nan,
                         'pool_size': pool_size,
                         'read_depth': read_depth,
                         'seq_error': seq_error,
                         'num_sites': int(1e6)},
                        ignore_index=True
                    )

                print('stat 2 pop')
                # statistical 2 pop
                for beta in np.linspace(0.1, 1., num=5):
                    true_freqs = np.zeros((int(1000), 2))
                    true_freqs_1 = np.random.beta(0.1, 0.5,
                                                  size=int(1000))
                    true_freqs_2 = np.random.beta(
                        beta * true_freqs_1 / (1 - true_freqs_1 + 1e-10),
                        beta
                    )
                    true_freqs[:, 0] = true_freqs_1
                    true_freqs[:, 1] = true_freqs_2
                    n0, n1 = statistical_model(true_freqs,
                                               pool_size,
                                               read_depth,
                                               seq_error)
                    fname = '{}_{}_{}_{}_stat_2_pop_beta_{}.txt'.format(
                        pool_size, read_depth, seq_error, rep, beta
                    )
                    print_count_matrix(n0, n1, fname)
                    res = res.append(
                        {'fname': fname,
                         'true_pairwise_het': np.nan,
                         'true_num_seg': np.nan,
                         'true_pi_w': pi_within(true_freqs),
                         'true_pi_b': pi_between(true_freqs),
                         'true_pi_t': pi_total(true_freqs),
                         'true_hudson_fst': fst_hudson(true_freqs),
                         'true_nei_fst': fst_nei(true_freqs),
                         'pool_size': pool_size,
                         'read_depth': read_depth,
                         'seq_error': seq_error,
                         'num_sites': int(1000)},
                        ignore_index=True
                    )

                print('constant size 1 pop')
                # constant size 1 pop
                for theta in np.linspace(1e-4, 1e-3, num=5):
                    (n0, n1), true_freqs, num_seg = msprime_model(
                        pool_size,
                        read_depth,
                        seq_error,
                        theta,
                        int(10e6),
                        1
                    )
                    fname = '{}_{}_{}_{}_const_1_pop_theta_{}.txt'.format(
                        pool_size, read_depth, seq_error, rep, theta
                    )
                    print_count_matrix(n0, n1, fname)
                    res = res.append(
                        {'fname': fname,
                         'true_pairwise_het': theta_pi(true_freqs),
                         'true_num_seg': num_seg,
                         'true_pi_w': np.nan,
                         'true_pi_b': np.nan,
                         'true_pi_t': np.nan,
                         'true_hudson_fst': np.nan,
                         'true_nei_fst': np.nan,
                         'pool_size': pool_size,
                         'read_depth': read_depth,
                         'seq_error': seq_error,
                         'num_sites': int(10e6)},
                        ignore_index=True
                    )

                print('exp growth 1 pop')
                # exponential 1 pop
                for theta in np.linspace(1e-4, 1e-3, num=5):
                    (n0, n1), true_freqs, num_seg = msprime_model(
                        pool_size,
                        read_depth,
                        seq_error,
                        theta,
                        int(10e6),
                        1,
                        exponential_growth=True
                    )
                    fname = '{}_{}_{}_{}_exp_growth_1_pop_theta_{}.txt'.format(
                        pool_size, read_depth, seq_error, rep, theta
                    )
                    print_count_matrix(n0, n1, fname)
                    res = res.append(
                        {'fname': fname,
                         'true_pairwise_het': theta_pi(true_freqs),
                         'true_num_seg': num_seg,
                         'true_pi_w': np.nan,
                         'true_pi_b': np.nan,
                         'true_pi_t': np.nan,
                         'true_hudson_fst': np.nan,
                         'true_nei_fst': np.nan,
                         'pool_size': pool_size,
                         'read_depth': read_depth,
                         'seq_error': seq_error,
                         'num_sites': int(10e6)},
                        ignore_index=True
                    )

                print('2 pop')
                # 2 pop for fst
                for t_div in np.linspace(0, 1, num=5):
                    (n0, n1), true_freqs = msprime_model(
                        pool_size,
                        read_depth,
                        seq_error,
                        0.001,
                        int(1e6),
                        2,
                        d_time=t_div
                    )
                    fname = '{}_{}_{}_{}_const_2_pop_t_div_{}.txt'.format(
                        pool_size, read_depth, seq_error, rep, t_div
                    )
                    print_count_matrix(n0, n1, fname)
                    res = res.append(
                        {'fname': fname,
                         'true_pairwise_het': np.nan,
                         'true_num_seg': np.nan,
                         'true_pi_w': pi_within(true_freqs),
                         'true_pi_b': pi_between(true_freqs),
                         'true_pi_t': pi_total(true_freqs),
                         'true_hudson_fst': fst_hudson(true_freqs),
                         'true_nei_fst': fst_nei(true_freqs),
                         'pool_size': pool_size,
                         'read_depth': read_depth,
                         'seq_error': seq_error,
                         'num_sites': int(1e6)},
                        ignore_index=True
                    )

res.to_csv('sim_metadata.csv', index=False)
