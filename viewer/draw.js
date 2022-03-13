// Trim whitespace from PDB data
// TODO: Do this as a pre-processing step.
molecule.atoms.forEach(a => {
  a.element = a.element.trim();
  a.residueName = a.residueName.trim();
  a.chain = a.chain.trim();
  a.name = a.name.trim();
})

// Properly sort atoms
molecule.atoms = molecule.atoms.sort((a, b) => a.serialNumber - b.serialNumber);

// Group by bases
const bases = (() => {
  let bases = {};
  molecule.atoms.forEach(a => {
    let seqN = a.residueSequenceNumber;
    if (bases[seqN] === undefined) {
      bases[seqN] = {
        type: a.residueName,
        residueSequenceNumber: a.residueSequenceNumber,
        chain: a.chain,
        atoms: {}
      };
    }
    // TODO: Check for dups?
    bases[seqN].atoms[a.name] = {x: a.x, y: a.y, z: a.z, element: a.element};
  })
  return Object.values(bases);
})();

const elementColors = {
  C: 0x808080, // gray
  N: 0x0000ff, // blue
  O: 0xff0000, // red
  P: 0xffa500, // orange
}

const baseColors = {
  A: 0xff0000, // red
  G: 0xff00ff, // fuchsia
  C: 0x00ff00, // lime
  U: 0xffff00, // yellow
}

// Atom is an object with cartesian coordinates
function sphereFromAtom(atom, r=0.5) {
  let g = new THREE.SphereGeometry(r, 20, 20);
  let m = new THREE.MeshToonMaterial();
  let s = new THREE.Mesh(g, m);
  s.position.set(atom.x, atom.y, atom.z);
  s.material.color.setHex(elementColors[atom.element]);
  return s;
}

// Phosphorus backbone
function drawBackbone(atoms, atm, color=0xffa500, w=5) {
  let points = [];
  molecule.atoms.forEach(a => {
    if (a.name === atm) {
      points.push(new THREE.Vector3(a.x, a.y, a.z));
    }
  })

  let m = new THREE.LineBasicMaterial({
    color: color, // orange by default
    linewidth: w
  });
  let g = new THREE.BufferGeometry().setFromPoints(points);
  let l = new THREE.Line(g, m);

  scene.add(l);
}

// Helper function which iterates through atoms in order
function _drawBonds(atoms, w, edgeList) {
  edgeList.forEach(edges => {
    let points = [];
    edges.forEach(e => {
      let a = atoms[e];
      points.push(new THREE.Vector3(a.x, a.y, a.z));
    })

    let m = new THREE.LineBasicMaterial( { color: 0xffffff, linewidth: w} );
    let g = new THREE.BufferGeometry().setFromPoints(points);
    let l = new THREE.Line(g, m);

    scene.add(l);
  })
}

function drawBonds(base, w) {
  _drawBonds(base.atoms, w, drawingOrder[base.type]);
}

// Draw nucleic acid bonds in proper order
const drawingOrder = {
  G: [
    ["C1'", "N9", "C8", "N7", "C5", "C6", "N1", "C2", "N3", "C4", "N9"],
    ["C4", "C5"],
    ["C6", "O6"],
    ["C2", "N2"],
    // Ribose
    ["P", "C5'", "C4'", "C3'", "C2'", "C1'"],
    ["C3'", "O3'"],
    ["C2'", "O2'"],
    ["C4'", "O4'", "C1'"]
  ],
  A: [
    ["C1'", "N9", "C8", "N7", "C5", "C6", "N1", "C2", "N3", "C4", "N9"],
    ["C4", "C5"],
    ["C6", "N6"],
    // Ribose
    ["P", "C5'", "C4'", "C3'", "C2'", "C1'"],
    ["C3'", "O3'"],
    ["C2'", "O2'"],
    ["C4'", "O4'", "C1'"]
  ],
  C: [
    ["C1'", "N1", "C6", "C5", "C4", "N3", "C2", "N1"],
    ["C4", "N4"],
    ["C2", "O2"],
    // Ribose
    ["P", "C5'", "C4'", "C3'", "C2'", "C1'"],
    ["C3'", "O3'"],
    ["C2'", "O2'"],
    ["C4'", "O4'", "C1'"]
  ],
  U: [
    ["C1'", "N1", "C6", "C5", "C4", "N3", "C2", "N1"],
    ["C4", "O4"],
    ["C2", "O2"],
    // Ribose
    ["P", "C5'", "C4'", "C3'", "C2'", "C1'"],
    ["C3'", "O3'"],
    ["C2'", "O2'"],
    ["C4'", "O4'", "C1'"]
  ]
}

// Draw all atoms as spheres
function drawSpheres(atoms, w=0.5) {
  atoms.forEach(a => {
    scene.add(sphereFromAtom(a, w));
  });
}

function drawHalfRung(base, w=5) {
  let t = base.type;
  let c3p = base.atoms["C3'"];

  var points = [new THREE.Vector3(c3p.x, c3p.y, c3p.z)];

  if (t === "A" || t === "G") { // Purines
    // C3' -> N1
    let n1 = base.atoms["N1"];
    points.push(new THREE.Vector3(n1.x, n1.y, n1.z));
  } else { // Pyrimidines
    // C3' -> N3
    let n3 = base.atoms["N3"];
    points.push(new THREE.Vector3(n3.x, n3.y, n3.z));
  }

  let m = new THREE.LineBasicMaterial({
    color: baseColors[t],
    linewidth: w
  });
  let g = new THREE.BufferGeometry().setFromPoints(points);
  let l = new THREE.Line(g, m);

  scene.add(l);
}

drawBackbone(molecule.atoms, "C3'");
bases.forEach(b => drawHalfRung(b));
// bases.forEach(b => drawBonds(b, 2));
// drawSpheres(molecule.atoms, 0.5);
