import { Component, OnInit, Input,ElementRef,ViewEncapsulation, ViewChild } from '@angular/core';
import * as d3 from 'd3';
import {BaseChartDirective} from 'ng2-charts/ng2-charts';
import {SimpleChanges} from '@angular/core';
export type Datum = {name: string, value: number};
@Component({
  encapsulation: ViewEncapsulation.None,
  selector: 'barcharts',
  templateUrl: './barcharts.component.html',
  styleUrls: ['./barcharts.component.css']
})

export class BarchartsComponent  {
@ViewChild(BaseChartDirective)
public chart: BaseChartDirective;
 public barChartOptions:any = {
    title: {
            display: true,
            text: 'Amino acid count by pocket'
        },
    scales: {
	yAxes: [{
	   scaleLabel: {
		display: true,
		labelString: 'Number of residues'
		}
	}],
	xAxes: [{
                display: true,
                labelString: 'Essential Amino Acids'
            }]
    }
  };
 
  
 
  public barChartLabels:string[] = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'];          
  public barChartType:string = 'bar';                                                                  
  public barChartLegend:boolean = true;                                                              
  
  public barChartData: any[];
  public sampleChartData:any[] = [
	{data: [65, 59, 80, 81, 56, 55, 40, 40, 30, 53, 30,
34, 65, 45, 35, 34, 89, 29 ,32, 39], label: 'Series A'},                                              
 	{data: [28, 48, 40, 19, 86, 27, 90, 80,
28, 49, 83, 92 ,73, 14, 90, 87, 17, 54, 38, 50], label: 'Series B'}
  ];                                                                                                 
 
  @Input() inputBarChartData : any[];
  

  ngOnChanges(){
	this.barChartData = this.inputBarChartData;
        this.redrawChart(); 
  }

 
  redrawChart(){
	 if (this.chart !== undefined) {
       		this.chart.chart.destroy();
       		this.chart.chart = 0;
       		this.chart.datasets = this.barChartData;
       		this.chart.labels = this.barChartLabels;
      		 this.chart.ngOnInit();
    	 }
   } 

  // events
  public chartClicked(e:any):void {
    console.log(e);
  }
 
  public chartHovered(e:any):void {
    console.log(e);
  }

 
  public randomize():void {
    // Only Change 3 values
    let data = [
      Math.round(Math.random() * 100),
      59,
      80,
      (Math.random() * 100),
      56,
      (Math.random() * 100),
      40];
    let clone = JSON.parse(JSON.stringify(this.barChartData));
    clone[0].data = data;
    this.barChartData = clone;
}
i}
