/* AUTO-GENERATED */
package JSci.maths.matrices;

import java.awt.Dimension;
import java.util.Hashtable;
import JSci.maths.algebras.Algebra;
import JSci.maths.fields.Ring;
import JSci.maths.groups.AbelianGroup;

public final class DoubleMatrixAlgebra implements Algebra, Ring {
	private static final Hashtable algebras = new Hashtable();
	static DoubleMatrixAlgebra get(int rows, int cols) {
		Dimension dim = new Dimension(rows, cols);
		DoubleMatrixAlgebra algebra = (DoubleMatrixAlgebra) algebras.get(dim);
		if(algebra == null) {
			algebra = new DoubleMatrixAlgebra(rows, cols);
			algebras.put(dim, algebra);
		}
		return algebra;
	}

	private final int rows;
	private final int cols;
	private AbstractDoubleMatrix zero;
	private AbstractDoubleSquareMatrix one;
	private DoubleMatrixAlgebra(int rows, int cols) {
		this.rows = rows;
		this.cols = cols;
	}
	/**
	 * Returns the (right) identity.
	 */
	public Ring.Member one() {
		if(one == null)
			one = DoubleDiagonalMatrix.identity(cols);
		return one;
	}
	public boolean isOne(Ring.Member r) {
		return one().equals(r);
	}
	public AbelianGroup.Member zero() {
		if(zero == null)
			zero = new DoubleMatrix(rows, cols);
		return zero;
	}
	public boolean isZero(AbelianGroup.Member r) {
		return zero().equals(r);
	}
	public boolean isNegative(AbelianGroup.Member a, AbelianGroup.Member b) {
		return zero().equals(a.add(b));
	}
}
