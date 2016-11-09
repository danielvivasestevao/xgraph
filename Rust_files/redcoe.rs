// Written using Rust v1.12.0

use std::collections::HashMap;
use std::collections::hash_map::Entry;
use std::collections::HashSet;
use std::fmt;
use std::fs::File;
use std::io::Read;
use std::io::Write;
use std::env;

/********************************* POINT *********************************/

#[derive(Debug, Clone, Copy)]
struct Point {
	x: f64,
	y: f64,
}

impl fmt::Display for Point {
	fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
		write!(f, "({}, {})", self.x, self.y)
	}
}

impl PartialEq for Point {
	fn eq(&self, other: &Point) -> bool {
		let x_comp: bool = &self.x == &other.x;
		let y_comp: bool = &self.y == &other.y;
		x_comp & y_comp
	}
}

/********************************* EDGE *********************************/

#[derive(Debug, Clone, Copy)]
struct Edge {
	left: u16,
	right: u16,
}

impl Edge {
	fn new(p1: u16, p2: u16, point_to_position: &HashMap<u16, Point>) -> Edge {
		let x1 = point_to_position.get(&p1).unwrap().x;
		let x2 = point_to_position.get(&p2).unwrap().x;
		let p_left: u16 = if x1 <= x2 {p1} else {p2};
		let p_right: u16 = if p_left == p1 {p2} else {p1};
		Edge { left: p_left, right: p_right}
	}
}

impl fmt::Display for Edge {
	fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
		write!(f, "({}, {})", self.left, self.right)
	}
}

impl PartialEq for Edge {
	fn eq(&self, other: &Edge) -> bool {
		let left_comp: bool = &self.left == &other.left;
		let right_comp: bool = &self.right == &other.right;
		left_comp & right_comp
	}
}

/******************************* TRIANGLE *******************************/

#[derive(Debug, Clone, Copy)]
struct Triangle {
	left: u16,
	middle: u16,
	right: u16
}

impl Triangle {
	fn new(p1: u16, p2: u16, p3: u16, point_to_position: &HashMap<u16, Point>) -> Triangle {
		let t1: (u16, f64) = (p1, point_to_position.get(&p1).unwrap().x);
		let t2: (u16, f64) = (p2, point_to_position.get(&p2).unwrap().x);
		let t3: (u16, f64) = (p3, point_to_position.get(&p3).unwrap().x);
		let mut v = [t1, t2, t3];
		v.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
		Triangle { left: v[0].0, middle: v[1].0, right: v[2].0 }
	}

	fn is_collinear(&self, point_to_position: &HashMap<u16, Point>) -> bool {
		let b1: bool = point_to_position.get(&self.left).unwrap().x == point_to_position.get(&self.middle).unwrap().x;
		let b2: bool = point_to_position.get(&self.middle).unwrap().x == point_to_position.get(&self.right).unwrap().x;

		let b3: bool = point_to_position.get(&self.left).unwrap().y == point_to_position.get(&self.middle).unwrap().y;
		let b4: bool = point_to_position.get(&self.middle).unwrap().y == point_to_position.get(&self.right).unwrap().y;

		return (b1 & b2) || (b3 & b4)
	}
}

impl fmt::Display for Triangle {
	fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
		write!(f, "({}, {}, {})", self.left, self.middle, self.right)
	}
}

impl PartialEq for Triangle {
	fn eq(&self, other: &Triangle) -> bool {
		let left_comp: bool = &self.left == &other.left;
		let mid_comp: bool = &self.middle == &other.middle;
		let right_comp: bool = &self.right == &other.right;
		left_comp & mid_comp & right_comp
	}
}

/********************************* MAIN *********************************/

fn main() {
	let args: Vec<_> = env::args().collect();
	// Do NOT use commas in the argument lists! Seperate all values with whitespace.
	// 1st arg: boolean which determines whether to perform redundancy or coexistence operations ("true"->red, "false"->coe)
	// 2nd arg: file with point ids and their positions as "id x_val y_val "
	// 3rd arg (if 1st arg is "true") : file with edges by node id; parse every edge as "id1 id2 "
	// 3rd arg (if 1st arg is "false"): file with list of triangles (triples of node ids) as "id1 id2 id3"
	// Example: >main.exe point2pos.txt edge_list_str.txt

	let red_or_coe: bool = args[1].parse().unwrap();

	let mut point_to_position_py = String::new();
	{
		let mut f = File::open(&args[2]).expect("Unable to open Point-to-Position file");
		f.read_to_string(&mut point_to_position_py).expect("Unable to read from Point-to-Position file");
	}

	// Getting the event list and id->pos HashMap
	let mut event_list: Vec<u16> = vec![];
	let mut point_to_position = HashMap::<u16, Point>::new();
	{
		let mut point_vector: Vec<(u16, Point)> = parse_tuples(&point_to_position_py);
		point_vector.sort_by(|a, b| a.1.x.partial_cmp(&b.1.x).unwrap());
		for tuple in point_vector.iter() {
			event_list.push(tuple.0);
			point_to_position.insert(tuple.0, tuple.1);
		}
	}

	// Only needed for intersection calc
	let mut edge_list_py = String::new();
	let mut edge_list: Vec<Edge> = vec![];
	if red_or_coe {
		let mut f = File::open(&args[3]).expect("Unable to open edge list file");
		f.read_to_string(&mut edge_list_py).expect("Unable to read from edge list file");
		edge_list = parse_edges(&edge_list_py, &point_to_position);
	}

	// Only needed for coexistence event calc
	let mut triangle_list_py = String::new();
	let mut triangles: Vec<Triangle> = vec![];
	if !red_or_coe {
		let mut f = File::open(&args[3]).expect("Unable to open triangle list file");
		f.read_to_string(&mut triangle_list_py).expect("Unable to read from triangle list file");
		triangles = parse_triangles(&triangle_list_py, &point_to_position);
		triangles.sort_by(
			|t1, t2|
			point_to_position.get(&t1.left).unwrap().x.partial_cmp(&point_to_position.get(&t2.right).unwrap().x).unwrap()
		);
	}

	if red_or_coe {
		let intersections: Vec<(Edge, Edge)> = simple_bentley_ottmann(event_list, edge_list, &point_to_position);
		println!("{:?}", intersections);
	} else {
		let coexistence_events: Vec<(u16, Triangle)> = coexistence_event_finder(event_list, triangles, &point_to_position);
		println!("{:?}", coexistence_events);
	}
}

/************************ SIMPLE BENTLEY-OTTMANN ************************/

fn simple_bentley_ottmann(event_list: Vec<u16>, edge_list: Vec<Edge>, point_to_position: &HashMap<u16, Point>) -> Vec<(Edge, Edge)> {
	// event_list must be SORTED by x values!

	let mut intersections: Vec<(Edge, Edge)> = vec![];

	// Getting HashMaps for point->{edges it starts, edges it ends}
	let mut point_to_start_edges = HashMap::<u16, Vec<Edge>>::new();
	let mut point_to_end_edges = HashMap::<u16, Vec<Edge>>::new();

	for edge in edge_list.iter() {
		let mut edge_vec = match point_to_start_edges.entry(edge.left) {
			Entry::Occupied(o) => o.into_mut(),
			Entry::Vacant(v) => v.insert(vec![]),
		};
		edge_vec.push(*edge);

		let mut edge_vec = match point_to_end_edges.entry(edge.right) {
			Entry::Occupied(o) => o.into_mut(),
			Entry::Vacant(v) => v.insert(vec![]),
		};
		edge_vec.push(*edge); 
		//point_to_end_edges.insert(edge.right, edge);
	}

	let mut active_edges: Vec<Edge> = vec![];
	//active_edges.push(Edge {left: 5, right: 9});
	//remove_edge(&mut active_edges, Edge {left: 5, right: 9});

	let mut c: u32 = 0;

	for point in event_list.iter() {
		let start_edge_empty: bool = match point_to_start_edges.entry(*point) {
			Entry::Vacant(_) => true,
			Entry::Occupied(_) => false
		};

		
		if !start_edge_empty {
			for edge1 in point_to_start_edges.get(&point).unwrap().iter() {
				for edge2 in active_edges.iter() {
					c = c + 1;
					if intersect(*edge1, *edge2, &point_to_position) {
						intersections.push((*edge1, *edge2));
					}
				}
				active_edges.push(*edge1);
			}
		}

		let end_edge_empty: bool = match point_to_end_edges.entry(*point) {
			Entry::Vacant(_) => true,
			Entry::Occupied(_) => false
		};

		if !end_edge_empty {
			for end_edge in point_to_end_edges.get(&point).unwrap().iter() {
				remove_edge(&mut active_edges, *end_edge);
			}
		}
	}

	// println!("Intersection checks: {}", c);
	intersections
}


fn remove_edge(vec: &mut Vec<Edge>, edge: Edge) {
	let index = vec.iter().position(|x| *x == edge).unwrap();
	vec.remove(index);
}


/************************** TRIANGLE-NODE EVENTS **************************/

fn coexistence_event_finder(node_list: Vec<u16>, triangle_list: Vec<Triangle>, point_to_position: &HashMap<u16, Point>) -> Vec<(u16, Triangle)> {
	// node_list must be SORTED by x values!

	let mut events: Vec<(u16, Triangle)> = vec![];

	for node in node_list.iter() {
		let node_pos: Point = *point_to_position.get(&node).unwrap();
		for triangle in triangle_list.iter().filter(
				|t|
				point_to_position.get(&t.right).unwrap().x > node_pos.x &&
				point_to_position.get(&t.left).unwrap().x < node_pos.x &&
				*node != t.left && *node != t.middle && *node != t.right
			) {
			// check if node lies within the triangle
			if is_point_in_triangle(node_pos, *point_to_position.get(&triangle.left).unwrap(), *point_to_position.get(&triangle.middle).unwrap(), *point_to_position.get(&triangle.right).unwrap()) {
				events.push((*node, *triangle));
			}
		}
	}

	events
}

/// Written with help from...
/// https://totologic.blogspot.de/2014/01/accurate-point-in-triangle-test.html
fn is_point_in_triangle(p: Point, t1: Point, t2: Point, t3: Point) -> bool {
	// if the point lies exactly on an edge of the triangle, it is disregarded
	if is_collinear(p, t1, t2) || is_collinear(p, t2, t3) || is_collinear(p, t1, t3) {
		return false;
	}

	let denominator: f64 = (t2.y - t3.y) * (t1.x - t3.x) + (t3.x - t2.x) * (t1.y - t3.y);

    /*
    The denominator is zero if all three nodes of the triangle lie on one line.
    In that case the area of the triangle is zero, and thus cannot have a node
    within itself.
    */

    if denominator == 0.0 {
    	return false;
    }

    let a: f64 = ((t2.y - t3.y) * (p.x - t3.x) + (t3.x - t2.x) * (p.y - t3.y)) / denominator;
    let b: f64 = ((t3.y - t1.y) * (p.x - t3.x) + (t1.x - t3.x) * (p.y - t3.y)) / denominator; 
    let c: f64 = 1.0 - a - b;

    let a1: bool = 0.0 <= a;
    let a2: bool = a <= 1.0;

    let b1: bool = 0.0 <= b;
    let b2: bool = b <= 1.0;

    let c1: bool = 0.0 <= c;
    let c2: bool = c <= 1.0;

    return (a1 & a2) && (b1 & b2) && (c1 & c2);

}

fn is_collinear(p1: Point, p2: Point, p3: Point) -> bool {
	let b1: bool = p1.x == p2.x;
	let b2: bool = p2.x == p3.x;

	let b3: bool = p1.y == p2.y;
	let b4: bool = p2.y == p3.y;

	return (b1 & b2) || (b3 & b4);
}

/*************************** VECTOR OPERATIONS ***************************/

/// Written with help from...
/// https://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect
/// https://github.com/pgkelley4/line-segments-intersect/blob/master/js/line-segments-intersect.js
/// Accessed at 22nd Oct 2016
fn intersect(edge1: Edge, edge2: Edge, point_to_position: &HashMap<u16, Point>) -> bool {
	let p1: Point = *point_to_position.get(&edge1.left).unwrap();
	let p2: Point = *point_to_position.get(&edge1.right).unwrap();
	let q1: Point = *point_to_position.get(&edge2.left).unwrap();
	let q2: Point = *point_to_position.get(&edge2.right).unwrap();

	// If the edges share any endpoints, they do not intersect
	if (edge1.left == edge2.left) || (edge1.left == edge2.right) ||
		(edge1.right == edge2.left) || (edge1.right == edge2.right) {
		return false;
	}
	
	let r: Point = subtract(p2, p1);
	let s: Point = subtract(q2, q1);

	let u_numerator: f64 = crossproduct(subtract(q1, p1), r);
	let denominator: f64 = crossproduct(r, s);

	if u_numerator == 0.0 && denominator == 0.0 {  // collinear
		// overlapping lines count as intersecting
		
		return !(( (q1.x - p1.x < 0.0) &&
				   (q1.x - p2.x < 0.0) &&
				   (q2.x - p1.x < 0.0) &&
				   (q2.x - p2.x < 0.0) ) ||
				 ( (q1.y - p1.y < 0.0) &&
				   (q1.y - p2.y < 0.0) &&
				   (q2.y - p1.y < 0.0) &&
				   (q2.y - p2.y < 0.0) ));
	}

	if denominator == 0.0 {
		// lines are parallel and non-intersecting
		return false;
	}

	// u = (p - q) x r / (r x s)
	let u: f64 = u_numerator / denominator;
	// t = (q - p) x s / (r x s)
	let t: f64 = crossproduct(subtract(*point_to_position.get(&edge2.left).unwrap(), *point_to_position.get(&edge1.left).unwrap()), s) / denominator;

	return (t >= 0.0) && (t <= 1.0) && (u >= 0.0) && (u <= 1.0);
}


fn subtract(point1: Point, point2: Point) -> Point {
	let x: f64 = point1.x - point2.x;
	let y: f64 = point1.y - point2.y;
	Point { x: x, y: y}
}


fn crossproduct(point1: Point, point2: Point) -> f64 {
	return point1.x * point2.y - point1.y * point2.x;
}

/***************************** WRITE TO FILE *****************************/

fn write_edge_vec(edge_list: Vec<(Edge, Edge)>) {
	let mut f = File::create("edge_list.txt").unwrap();
	for tup in edge_list.iter() {
		f.write("(".as_bytes()).unwrap();
		f.write(tup.0.to_string().as_bytes()).unwrap();
		f.write(",".as_bytes()).unwrap();
		f.write(tup.1.to_string().as_bytes()).unwrap();
		f.write(")\n".as_bytes()).unwrap();
	}
}

/******************************** PARSING ********************************/

fn parse_tuples(arg_str: &String) -> Vec<(u16, Point)> {
	let mut result: Vec<(u16, Point)> = vec![];
	let mut arg_iter = arg_str.split_whitespace();

	loop {
		let id: u16 = match arg_iter.next() {
			None => break, 
			Some(v) => v.parse::<u16>().expect("parse error on id"),
		};
		let x_val: f64 = arg_iter.next().unwrap().parse::<f64>().expect("parse error on first float");
		let y_val: f64 = arg_iter.next().unwrap().parse::<f64>().expect("parse error on second float");
		result.push((id, Point { x: x_val, y: y_val }));
	}

	result
}

fn parse_edges(arg_str: &String, point_to_position: &HashMap<u16, Point>) -> Vec<Edge> {
	let mut result: Vec<Edge> = vec![];
	let mut arg_iter = arg_str.split_whitespace();

	loop {
		let p1: u16 = match arg_iter.next() {
			None => break,
			Some(v) => v.parse::<u16>().expect("parse error on first point of edge"),
		};
		let p2: u16 = arg_iter.next().unwrap().parse::<u16>().expect("parse error on second point of edge");
		result.push(Edge::new(p1, p2, &point_to_position));
	}

	result
}

fn parse_triangles(arg_str: &String, point_to_position: &HashMap<u16, Point>) -> Vec<Triangle> {
	let mut result: Vec<Triangle> = vec![];
	let mut arg_iter = arg_str.split_whitespace();

	loop {
		let p1: u16 = match arg_iter.next() {
			None => break,
			Some(v) => v.parse::<u16>().expect("parse error on first point of triangle"),
		};
		let p2: u16 = arg_iter.next().unwrap().parse::<u16>().expect("parse error on second point of triangle");
		let p3: u16 = arg_iter.next().unwrap().parse::<u16>().expect("parse error on third point of triangle");
		result.push(Triangle::new(p1, p2, p3, &point_to_position));
	}

	result	
}