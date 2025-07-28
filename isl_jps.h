#ifndef ISL_JPS_H_
#define ISL_JPS_H_

/* isl_jps - v0.1 

   Public domain grid pathfinder no warranty implied; use at your own risk.

   Do this:
       #define ISL_JPS_IMPLEMENTATION
   before you include this file in *one* C or C++ file to create the
   implementation.

   To static link also add:
       #define ISLJ_STATIC

	 QUICK NOTES:

   This is simple wrapper around library written by Ari Rahikkala, taken from
   GitHub: https://github.com/arirahikkala/astar-jps, original README with
   LICENSE is at the of this file.
   
   LICENSE:
     See end of file for license information.
*/

#ifndef ISLJ_DEF
#ifdef ISLJ_STATIC
#define ISLJ_DEF static
#else
#define ISLJ_DEF extern
#endif
#endif

/* Run A* pathfinding over uniform-cost 2d grids using jump point search.

   grid: 0 if obstructed, non-0 if non-obstructed(the value is ignored beyond that).
   solLength: out-parameter, will contain the solution length
   boundX: width of the grid
   boundY: height of the grid
   start: index of the starting node in the grid
   end: index of the ending node in the grid

   return value: Array of node indexes making up the solution, of length solLength, in reverse order
 */

#ifdef __cplusplus
extern "C" {
#endif

ISLJ_DEF int *islj_compute(const char *grid, int *solLength, int boundX, int boundY, int start, int end);
ISLJ_DEF int *islj_unopt_compute(const char *grid, int *solLength, int boundX, int boundY, int start, int end);
ISLJ_DEF int islj_index(int width, int x, int y); /* Compute cell indexes from cell coordinates and the grid width */
ISLJ_DEF void islj_coord(int width, int node, int *x, int *y);  /* Compute coordinates from a cell index and the grid width */

#ifdef __cplusplus
}
#endif
#endif // ISL_JPS_H_

#ifdef ISL_JPS_IMPLEMENTATION

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

typedef struct coord {
	int x;
	int y;
} islj_xy;

typedef struct islj_pqueue_item {
	double priority;
	int value;
} islj_pqueue_item;

typedef struct islj_pqueue {
	int size;
	unsigned int allocated;
	islj_pqueue_item *root;
	int *index;
	unsigned int indexAllocated;
} islj_pqueue;

static void islj_pq_insert(islj_pqueue *q, int value, double priority);
static void islj_pq_delete_min(islj_pqueue *q);
static islj_pqueue_item *islj_pq_find_min(const islj_pqueue *q);
static void islj_pq_change_priority(islj_pqueue *q, int ind, double newPriority);
static int islj_pq_priority_of(const islj_pqueue *q, int ind);
static int islj_pq_exists(const islj_pqueue *q, int ind);
static islj_pqueue *islj_pq_create(void);
static void islj_pq_free(islj_pqueue *q);
static int islj_pq_make_space(islj_pqueue *q, size_t size);
static int smallest_power_of_2_after(int);
static int islj_pq_place_at_end(islj_pqueue *q, islj_pqueue_item islj_pqueue_item);
static void islj_pq_sift_up(islj_pqueue *q, int i);
static void islj_pq_sift_down(islj_pqueue *q, int i);

typedef int node;

typedef struct astar {
	const char *grid;
	islj_xy bounds;
	node start;
	node goal;
	islj_pqueue *open;
	char *closed;
	double *gScores;
	node *cameFrom;
	int *solutionLength;
} islj_astar;

// The order of islj_directions is: 
// N, NE, E, SE, S, SW, W, NW 
typedef unsigned char islj_direction;
#define NO_DIRECTION 8
typedef unsigned char islj_directionset;

static double islj_estimate_distance(islj_xy start, islj_xy end);
static double islj_precise_distance(islj_xy start, islj_xy end);

static int islj_astar_init(islj_astar* astar, const char *grid, int *solLength, int boundX, int boundY, int start, int end);
static int islj_astar_is_enterable(islj_astar *astar, islj_xy coord);
static void islj_astar_add_to_open_set(islj_astar *astar, int node,  int nodeFrom);
static islj_directionset islj_astar_forced_neighbours(islj_astar *astar, islj_xy coord, islj_direction dir);
static islj_direction islj_astar_direction_we_came_from(islj_astar *astar, int node, int nodeFrom);
static int *islj_astar_record_solution(islj_astar *astar);
static int islj_astar_next_node_in_solution(islj_astar *astar, int *target, int node);
static int islj_astar_jump(islj_astar *astar, islj_direction dir, int start);

static islj_direction islj_next_direction_in_set(islj_directionset *dirs);
static islj_directionset islj_add_direction_to_set(islj_directionset dirs, islj_direction dir);
static islj_direction islj_direction_of_move(islj_xy to, islj_xy from);
static int islj_get_index(islj_xy bounds, islj_xy c);
static islj_xy islj_get_coord(islj_xy bounds, int c);
static int islj_contained(islj_xy bounds, islj_xy c);
static int islj_direction_is_diagonal(islj_direction dir);
static islj_xy islj_adjust_in_direction(islj_xy c, int dir);
static int implies(int a, int b);
static islj_directionset islj_natural_neighbours(islj_direction dir);

int smallest_power_of_2_after(int x) {
	if(x < 0) return 0;
	x |= x >> 1;
	x |= x >> 2;
	x |= x >> 4;
	x |= x >> 8;
	x |= x >> 16;
	return x+1;
}

int islj_pq_make_space(islj_pqueue *q, size_t size) {
	unsigned int newAllocated = smallest_power_of_2_after(size * sizeof(islj_pqueue_item));

	if(newAllocated <= q->allocated)
		return 0;

	q->root = realloc(q->root, newAllocated);
	if(NULL == q->root)
		exit(1); // I'm a bloody Haskell programmer, I don't have any business trying to do anything fancy when realloc fails

	q->allocated = newAllocated;
	return newAllocated;
}

int islj_pq_place_at_end(islj_pqueue *q, islj_pqueue_item islj_pqueue_item) {
	islj_pq_make_space(q, q->size + 1);
	q->root[q->size] = islj_pqueue_item;
	return q->size++;
}

void islj_pq_sift_up(islj_pqueue *q, int i) {
	if(0 == i)
		return;

	int p =(i - 1) / 2;

	if(q->root[p].priority < q->root[i].priority)
		return;

	q->index[q->root[i].value] = p;
	q->index[q->root[p].value] = i;

	islj_pqueue_item swap = q->root[i];
	q->root[i] = q->root[p];
	q->root[p] = swap;

	return islj_pq_sift_up(q, p);
}

void islj_pq_insert(islj_pqueue *q, int value, double pri) {
	islj_pqueue_item i;
	i.value = value;
	i.priority = pri;

	int newAllocated = smallest_power_of_2_after((value + 1) * sizeof(int));

	if((value + 1) * sizeof(int) > q->indexAllocated) {
		q->index = realloc(q->index, newAllocated);
		if(NULL == q->index)
			exit(1);
		for(unsigned int j = q->indexAllocated / sizeof(int); j < newAllocated / sizeof(int); j++)
			q->index[j] = -1;
		q->indexAllocated = newAllocated;
	}

	int p = islj_pq_place_at_end(q, i);

	q->index[q->root[p].value] = p;	

	islj_pq_sift_up(q, p);
}

void islj_pq_sift_down(islj_pqueue *q, int i) {
	int c = 1 + 2 * i;
	if(c >= q->size)
		return;

	if((c + 1 < q->size) &&(q->root[c].priority > q->root[(c+1)].priority))
		c++;

	if(q->root[i].priority < q->root[c].priority)
		return;

	q->index[q->root[c].value] = i;
	q->index[q->root[i].value] = c;

	islj_pqueue_item swap = q->root[i];
	q->root[i] = q->root[c];
	q->root[c] = swap;

	return islj_pq_sift_down(q, c);
}

void islj_pq_delete_min(islj_pqueue *q) {
	if(0 == q->size)
		return;

	q->index[q->root[0].value] = -1;
	q->size--;

	if(0 == q->size)
		return;

	q->index[q->root[q->size].value] = 0;
	q->root[0] = q->root[q->size];

	islj_pq_sift_down(q, 0);
} 

islj_pqueue_item *islj_pq_find_min(const islj_pqueue *q) {
	return q->root;
}

void islj_pq_change_priority(islj_pqueue *q, int ind, double newPriority) {
	int oldPriority = q->root[q->index[ind]].priority;
	q->root[q->index[ind]].priority = newPriority;
	if(oldPriority < newPriority)
		islj_pq_sift_down(q, q->index[ind]);
	else if(oldPriority > newPriority)
		islj_pq_sift_up(q, q->index[ind]);
}

void delete(islj_pqueue *q, int ind) {
	islj_pq_change_priority(q, ind, INT_MIN);
	islj_pq_delete_min(q);
}

int islj_pq_priority_of(const islj_pqueue *q, int ind) {
	return q->root[q->index[ind]].priority;
}

int islj_pq_exists(const islj_pqueue *q, int ind) {
	return (q->indexAllocated / sizeof(int) > ind) &&
		(-1 != q->index[ind]) && 
		(q->size > q->index[ind]);
}

islj_pqueue* islj_pq_create(void) {
	islj_pqueue *rv =(islj_pqueue*) malloc(sizeof(islj_pqueue));
	if(NULL == rv)
		exit(1);
	rv->size = 0;
	rv->allocated = 0;
	rv->root = NULL;
	rv->index = NULL;
	rv->indexAllocated = 0;
	return rv;
}

void islj_pq_free(islj_pqueue* q) {
	free(q->root);
	free(q->index);
	free(q);
}

// Distance metrics, you might want to change these to match your game mechanics

// Chebyshev distance metric for distance estimation by default
double islj_estimate_distance(islj_xy start, islj_xy end) {
	return fmax(abs(start.x - end.x), abs(start.y - end.y));
}

// Since we only work on uniform-cost maps, this function only needs
// to see the coordinates, not the map itself.
// Euclidean geometry by default.
// Note that since we jump over points, we actually have to compute 
// the entire distance - despite the uniform cost we can't just collapse
// all costs to 1
double islj_precise_distance(islj_xy start, islj_xy end) {
	if(start.x - end.x != 0 && start.y - end.y != 0)
		return sqrt(pow(start.x - end.x, 2) + 
				pow(start.y - end.y, 2)) ;
	else
		return abs(start.x - end.x) + abs(start.y - end.y);
}

// Below this point, not a lot that there should be much need to change!


// return and remove a islj_direction from the set
// returns NO_DIRECTION if the set was empty
islj_direction islj_next_direction_in_set(islj_directionset *dirs) {
	for(int i = 0; i < 8; i++) {
		char bit = 1 << i;
		if(*dirs & bit) {
			*dirs ^= bit;
			return i;
		}
	}
	return NO_DIRECTION;
}

islj_directionset islj_add_direction_to_set(islj_directionset dirs, islj_direction dir) {
	return dirs | 1 << dir;
}

/* Coordinates are represented either as pairs of an x-coordinate and
	 y-coordinate, or map indexes, as appropriate. islj_get_index and islj_get_coord
	 convert between the representations. */
int islj_get_index(islj_xy bounds, islj_xy c) {
	return c.x + c.y * bounds.x;
}

islj_xy islj_get_coord(islj_xy bounds, int c) {
	islj_xy rv = { c % bounds.x, c / bounds.x };
	return rv;
}

// is this coordinate islj_contained within the map bounds?
int islj_contained(islj_xy bounds, islj_xy c) {
	return c.x >= 0 && c.y >= 0 && c.x < bounds.x && c.y < bounds.y;
}

// is this coordinate within the map bounds, and also walkable?
int islj_astar_is_enterable(islj_astar *astar, islj_xy coord) {
	node node = islj_get_index(astar->bounds, coord);
	return islj_contained(astar->bounds, coord) && astar->grid[node];
}

int islj_direction_is_diagonal(islj_direction dir) {
	return(dir % 2) != 0;
}

// the coordinate one tile in the given islj_direction
islj_xy islj_adjust_in_direction(islj_xy c, int dir) {
	// we want to implement "rotation" - that is, for instance, we can
	// subtract 2 from the islj_direction "north" and get "east"
	// C's modulo operator doesn't quite behave the right way to do this,
	// but for our purposes this kluge should be good enough
	switch((dir + 65536) % 8) {
		case 0: return(islj_xy) {c.x, c.y - 1};
		case 1: return(islj_xy) {c.x + 1, c.y - 1};
		case 2: return(islj_xy) {c.x + 1, c.y };
		case 3: return(islj_xy) {c.x + 1, c.y + 1};
		case 4: return(islj_xy) {c.x, c.y + 1};
		case 5: return(islj_xy) {c.x - 1, c.y + 1};
		case 6: return(islj_xy) {c.x - 1, c.y};
		case 7: return(islj_xy) {c.x - 1, c.y - 1};
	}
	return(islj_xy) { -1, -1 };
}

// logical implication operator
int implies(int a, int b) {
	return a ? b : 1;	
}

/* Harabor's explanation of exactly how to determine when a cell has forced
	 neighbours is a bit unclear IMO, but this is the best explanation I could
	 figure out. I won't go through everything in the paper, just the extra
	 insights above what I thought was immediately understandable that it took
	 to actually implement this function.

	 First, to introduce the problem, we're looking at the immedate neighbours
	 of a cell on the grid, considering what tile we arrived from.

	 ...  This is the basic situation we're looking at. Supposing the top left
	 -X.  period is cell(0,0), we're coming in to cell(1, 1) from(0, 1).
	 ...  

	 ...  The other basic case, the diagonal case. All other cases are obviously
	 .X.  derivable from these two cases by symmetry.
	 /..

	 The question is: Given that some tiles might have walls, *how many tiles
	 are there that we can reach better by going through the center tile than
	 any other way?*(for the horizontal case it's ok to only be able to reach
	 them as well some other as through the center tile too)

	 In case there are no obstructions, the answers are simple: In the horizontal
	 or vertical case, the cell directly ahead; in the diagonal case, the three
	 cells ahead.

	 The paper is pretty abstract about what happens when there *are* 
	 obstructions, but fortunately the abstraction seems to collapse into some
	 fairly simple practical cases:

	 123  Position 4 is a natural neighbour(according to the paper's terminology)
	 -X4  so we don't need to count it. Positions 1, 2, 5 and 6 are accessible
	 567  without going through the center tile. This leaves positions 3 and 7
	 to be looked at.

	 Considering position 3(everything here also follows for 7 by symmetry):
	 If 3 is obstructed, then it doesn't matter what's in position in 2.
	 If 3 is free and 2 is obstructed, 3 is a forced neighbour.
	 If 3 is free and 2 is free, 3 is pruned(not a forced neighbour)

	 i.e. logically, 
	 3 is not a forced neighbour iff(3 is obstructed) implies(2 is obstructed).

	 Similar reasoning applies for the diagonal case, except with bigger angles.

*/
/*
	 static int hasForcedNeighbours(islj_astar *astar, islj_xy coord, int dir)
	 {
#define ENTERABLE(n) islj_astar_is_enterable(astar, \
islj_adjust_in_direction(coord, dir +(n)))
if(islj_direction_is_diagonal(dir))
return !implies(ENTERABLE(-2), ENTERABLE(-3)) ||
!implies(ENTERABLE(2), ENTERABLE(3));
else 
return !implies(ENTERABLE(-1), ENTERABLE(-2)) ||
!implies(ENTERABLE(1), ENTERABLE(2));
#undef ENTERABLE
}
*/
islj_directionset islj_astar_forced_neighbours(islj_astar *astar, islj_xy coord, islj_direction dir) {
	if(dir == NO_DIRECTION)
		return 0;

	islj_directionset dirs = 0;
#define ENTERABLE(n) islj_astar_is_enterable(astar, \
		islj_adjust_in_direction(coord,(dir +(n)) % 8))
	if(islj_direction_is_diagonal(dir)) {
		if(!implies(ENTERABLE(6), ENTERABLE(5)))
			dirs = islj_add_direction_to_set(dirs,(dir + 6) % 8);
		if(!implies(ENTERABLE(2), ENTERABLE(3)))
			dirs = islj_add_direction_to_set(dirs,(dir + 2) % 8);
	}
	else {
		if(!implies(ENTERABLE(7), ENTERABLE(6)))
			dirs = islj_add_direction_to_set(dirs,(dir + 7) % 8);
		if(!implies(ENTERABLE(1), ENTERABLE(2)))
			dirs = islj_add_direction_to_set(dirs,(dir + 1) % 8);
	}	
#undef ENTERABLE	
	return dirs;
}

islj_directionset islj_natural_neighbours(islj_direction dir) {
	if(dir == NO_DIRECTION)
		return 255;

	islj_directionset dirs = 0;
	dirs = islj_add_direction_to_set(dirs, dir);
	if(islj_direction_is_diagonal(dir)) {
		dirs = islj_add_direction_to_set(dirs,(dir + 1) % 8);
		dirs = islj_add_direction_to_set(dirs,(dir + 7) % 8);
	}

	return dirs;
}

void islj_astar_add_to_open_set(islj_astar *astar, int node,  int nodeFrom) {
	islj_xy nodeCoord = islj_get_coord(astar->bounds, node);
	islj_xy nodeFromCoord = islj_get_coord(astar->bounds, nodeFrom);

	if(!islj_pq_exists(astar->open, node)) {
		astar->cameFrom[node] = nodeFrom;
		astar->gScores[node] = astar->gScores[nodeFrom] + islj_precise_distance(nodeFromCoord, nodeCoord);
		islj_pq_insert(astar->open, node, astar->gScores[node] + islj_estimate_distance(nodeCoord, islj_get_coord(astar->bounds, astar->goal)));
	}
	else if(astar->gScores[node] > astar->gScores[nodeFrom] + islj_precise_distance(nodeFromCoord, nodeCoord)) {
		astar->cameFrom[node] = nodeFrom;
		int oldGScore = astar->gScores[node];
		astar->gScores[node] = astar->gScores[nodeFrom] + islj_precise_distance(nodeFromCoord, nodeCoord);
		double newPri = islj_pq_priority_of(astar->open, node) - oldGScore + astar->gScores[node];
		islj_pq_change_priority(astar->open, node, newPri);
	}	
}

// directly translated from "algorithm 2" in the paper
int islj_astar_jump(islj_astar *astar, islj_direction dir, int start) {
	islj_xy coord = islj_adjust_in_direction(islj_get_coord(astar->bounds, start), dir);
	int node = islj_get_index(astar->bounds, coord);
	if(!islj_astar_is_enterable(astar, coord))
		return -1;

	if(node == astar->goal || 
			islj_astar_forced_neighbours(astar, coord, dir)) {
		return node;
	}

	if(islj_direction_is_diagonal(dir)) {
		int next = islj_astar_jump(astar,(dir + 7) % 8, node);
		if(next >= 0)
			return node;

		next = islj_astar_jump(astar,(dir + 1) % 8, node);
		if(next >= 0)
			return node;
	}

	return islj_astar_jump(astar, dir, node);
}

// path interpolation between jump points in here
int islj_astar_next_node_in_solution(islj_astar *astar, int *target, int node) {
	islj_xy c = islj_get_coord(astar->bounds, node);
	islj_xy cTarget = islj_get_coord(astar->bounds, *target);

	if(c.x < cTarget.x) 
		c.x++;
	else if(c.x > cTarget.x)
		c.x--;

	if(c.y < cTarget.y) 
		c.y++;
	else if(c.y > cTarget.y)
		c.y--;

	node = islj_get_index(astar->bounds, c);

	if(node == *target)
		*target = astar->cameFrom[*target];

	return node;
}

// a bit more complex than the usual A* solution-recording method,
// due to the need to interpolate path chunks
int *islj_astar_record_solution(islj_astar *astar) {
	int rvLen = 1;
	*astar->solutionLength = 0;
	int target = astar->goal;
	int *rv = malloc(rvLen * sizeof(int));
	int i = astar->goal;

	for(;;) {
		i = islj_astar_next_node_in_solution(astar, &target, i);
		rv[*astar->solutionLength] = i;
		(*astar->solutionLength)++;
		if(*astar->solutionLength >= rvLen) {
			rvLen *= 2;
			rv = realloc(rv, rvLen * sizeof(int));
			if(!rv)
				return NULL;
		}
		if(i == astar->start)
			break;
	}

	(*astar->solutionLength)--; // don't include the starting tile
	return rv;
}


static islj_direction islj_direction_of_move(islj_xy to, islj_xy from) {
	if(from.x == to.x) {
		if(from.y == to.y)
			return -1;
		else if(from.y < to.y)
			return 4;
		else // from.y > to.y
			return 0;
	}
	else if(from.x < to.x) {
		if(from.y == to.y)
			return 2;
		else if(from.y < to.y)
			return 3;
		else // from.y > to.y
			return 1;
	}
	else { // from.x > to.x
		if(from.y == to.y)
			return 6;
		else if(from.y < to.y)
			return 5;
		else // from.y > to.y
			return 7;
	}
}

islj_direction islj_astar_direction_we_came_from(islj_astar *astar, int node, int nodeFrom) {
	if(nodeFrom == -1)
		return NO_DIRECTION;

	return islj_direction_of_move(islj_get_coord(astar->bounds, node), islj_get_coord(astar->bounds, nodeFrom));
}

int islj_astar_init(islj_astar* astar, const char *grid, int *solLength, int boundX, int boundY, int start, int end) {
	*solLength = -1;
	islj_xy bounds = {boundX, boundY};

	int size = bounds.x * bounds.y;

	if(start >= size || start < 0 || end >= size || end < 0)
		return 0;

	islj_xy startCoord = islj_get_coord(bounds, start);
	islj_xy endCoord = islj_get_coord(bounds, end);

	if(!islj_contained(bounds, startCoord) || !islj_contained(bounds, endCoord))
		return 0;

	astar->solutionLength = solLength;
	astar->bounds = bounds;
	astar->start = start;
	astar->goal = end;
	astar->grid = grid;

	astar->open = islj_pq_create();
	if(!astar->open)
		return 0;

	astar->closed = malloc(size);
	if(!astar->closed) {
		islj_pq_free(astar->open);
		return 0;
	}

	astar->gScores = malloc(size * sizeof(double));
	if(!astar->gScores) {
		islj_pq_free(astar->open);
		free(astar->closed);
		return 0;
	}

	astar->cameFrom = malloc(size * sizeof(int));
	if(!astar->cameFrom) {
		islj_pq_free(astar->open);
		free(astar->closed);
		free(astar->gScores);
		return 0;
	}

	memset(astar->closed, 0, size);

	astar->gScores[start] = 0;
	astar->cameFrom[start] = -1;

	islj_pq_insert(astar->open, astar->start, islj_estimate_distance(startCoord, endCoord));

	return 1;
}

int *islj_compute(const char *grid, int *solLength, int boundX, int boundY, int start, int end) {
	islj_astar astar;
	if(!islj_astar_init(&astar, grid, solLength, boundX, boundY, start, end))
		return NULL;

	islj_xy bounds = {boundX, boundY};
	islj_xy endCoord = islj_get_coord(bounds, end);

	while(astar.open->size) {
		int node = islj_pq_find_min(astar.open)->value; 
		islj_xy nodeCoord = islj_get_coord(bounds, node);
		if(nodeCoord.x == endCoord.x && nodeCoord.y == endCoord.y) {
			islj_pq_free(astar.open);
			free(astar.closed);
			free(astar.gScores);

			int *rv = islj_astar_record_solution(&astar);

			free(astar.cameFrom);

			return rv;
		}

		islj_pq_delete_min(astar.open);
		astar.closed[node] = 1;

		islj_direction from = islj_astar_direction_we_came_from(&astar, node, astar.cameFrom[node]);
		islj_directionset dirs = islj_astar_forced_neighbours(&astar, nodeCoord, from) | islj_natural_neighbours(from);

		for(int dir = islj_next_direction_in_set(&dirs); dir != NO_DIRECTION; dir = islj_next_direction_in_set(&dirs)) {
			int newNode = islj_astar_jump(&astar, dir, node);
			islj_xy newCoord = islj_get_coord(bounds, newNode);

			// this'll also bail out if jump() returned -1
			if(!islj_contained(bounds, newCoord))
				continue;

			if(astar.closed[newNode])
				continue;

			islj_astar_add_to_open_set(&astar, newNode, node);
		}
	}
	islj_pq_free(astar.open);
	free(astar.closed);
	free(astar.gScores);
	free(astar.cameFrom);

	return NULL;
}

int *islj_unopt_compute(const char *grid, int *solLength, int boundX, int boundY, int start, int end) {
	islj_astar astar;

	if(!islj_astar_init(&astar, grid, solLength, boundX, boundY, start, end))
		return NULL;

	islj_xy bounds = {boundX, boundY};
	islj_xy endCoord = islj_get_coord(bounds, end);

	while(astar.open->size) {
		int node = islj_pq_find_min(astar.open)->value; 
		islj_xy nodeCoord = islj_get_coord(bounds, node);
		if(nodeCoord.x == endCoord.x && nodeCoord.y == endCoord.y) {
			islj_pq_free(astar.open);
			free(astar.closed);
			free(astar.gScores);

			int *rv = islj_astar_record_solution(&astar);
			free(astar.cameFrom);

			return rv;
		}

		islj_pq_delete_min(astar.open);
		astar.closed[node] = 1;

		for(int dir = 0; dir < 8; dir++)
		{
			islj_xy newCoord = islj_adjust_in_direction(nodeCoord, dir);
			int newNode = islj_get_index(bounds, newCoord);

			if(!islj_contained(bounds, newCoord) || !grid[newNode])
				continue;

			if(astar.closed[newNode])
				continue;

			islj_astar_add_to_open_set(&astar, newNode, node);
		}
	}
	islj_pq_free(astar.open);
	free(astar.closed);
	free(astar.gScores);
	free(astar.cameFrom);
	return NULL;
}

int islj_index(int width, int x, int y) {
	return x + y * width;
}

void islj_coord(int width, int node, int *x, int *y) {
	*x = node % width;
	*y = node / width;
}
/*
------------------------------------------------------------------------------
This software is available under 2 licenses -- choose whichever you prefer.
------------------------------------------------------------------------------
ALTERNATIVE A - MIT License
Copyright (c) 2025 Ilya Kolbin
Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
------------------------------------------------------------------------------
ALTERNATIVE B - Public Domain (www.unlicense.org)
This is free and unencumbered software released into the public domain.
Anyone is free to copy, modify, publish, use, compile, sell, or distribute this
software, either in source code form or as a compiled binary, for any purpose,
commercial or non-commercial, and by any means.
In jurisdictions that recognize copyright laws, the author or authors of this
software dedicate any and all copyright interest in the software to the public
domain. We make this dedication for the benefit of the public at large and to
the detriment of our heirs and successors. We intend this dedication to be an
overt act of relinquishment in perpetuity of all present and future rights to
this software under copyright law.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
------------------------------------------------------------------------------
*/

/* Original README with LICENSE https://github.com/iskolbin/isl_astar/blob/master/README.md 
A simple C library for A* pathfinding over uniform-cost 2-dimensional grids. Meant to be embedded and modified as needed. See AStar.h for pertinent documentation, and TestAStar.c for a simple usage example.

Based on the well-known A* and binary heap algorithms, with jump point search from D. Harabor and A. Grastien. Online Graph Pruning for Pathfinding on Grid Maps. In National Conference on Artificial Intelligence (AAAI), 2011. Or, for those who of us who prefer clicking on links to tracking down academical references: http://grastien.net/ban/articles/hg-aaai11.pdf

Copyright 2011 Ari Rahikkala. All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are
permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice, this list of
      conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice, this list
      of conditions and the following disclaimer in the documentation and/or other materials
      provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE AUTHOR ''AS IS'' AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#endif // ISL_JPS_IMPLEMENTATION

