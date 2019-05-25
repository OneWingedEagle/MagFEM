package components;


import java.awt.Dimension;
import java.awt.Font;
import javax.swing.JTextField;


public class TextField extends JTextField {
	
	int size=(int)(10.0*GUI.screenWidth/1200);
	Font tfFont = new Font("Arial", 0, size);
	
	public TextField(){
	super();	
	
	setFont(tfFont);
	}
	public TextField(int size ){
		super();	
		setSize(new Dimension(10,size));
		setFont(tfFont);
		}
	public TextField(String str ){
		super();	
		setText(str);
		setFont(tfFont);
		}
	public TextField(int size, String str ){
		super();	
		setText(str);
		setSize(new Dimension(10,size));
		setFont(tfFont);
		}
}
